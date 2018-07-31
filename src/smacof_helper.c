/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file smacof_helper.c
 * @brief Implementiert die Funktionsprototypen aus smacof_helper.h der für den
 * SMACOF Algorithmus benötigten Funktionen.
 *
 * Diese Datei enthält die Implementationen der Funktionsprototypen der
 * benötigten Funktionen, zur Formulierung des SMACOF Algorithmus.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#include "smacof_helper.h"
#include "cuda_disfuncs.cuh"
#include "cuda_linalg.cuh"
#include "cuda_smacof_helper.cuh"
#include "disfuncs.h"
#include "linalg.h"
#include "m.h"
#include "tm.h"
#include "utils.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

m_t *guttmanTransformation(const tm_t *delta, const m_t *z, const tm_t *w,
                           const tm_t *pinv) {
  m_t *update;
  const real_t norm = 1.0 / delta->size;
  // OPENMP C version of guttman-transformation
  tm_t *dz = calcDeltaMatrix(z, EUCLIDEAN, 2);
  tm_t *bz = computeBofZ(delta, dz, w);
  m_t *bztimesz = TMmultM(bz, z);
  if (w == NULL) {
    // formula in 'linear' math: X_update = n^-1 * B(Z) * Z
    update = Mmultscalar(bztimesz, norm);
  } else {
    update = TMmultM(pinv, bztimesz);
  }
  // cleanup temporary stuff
  freeTM(dz);
  dz = NULL;
  freeTM(bz);
  bz = NULL;
  freeM(bztimesz);
  bztimesz = NULL;
  return update;
}

m_t *guttmanTransformation_nomem(const tm_t *delta, m_t *z, m_t *x,
                                 const tm_t *w, const tm_t *pinv,
                                 const tm_t *dz, tm_t *bz) {
  const real_t norm = 1.0 / delta->size;
  // OPENMP C version of guttman-transformation
  computeBofZ_nomem(delta, dz, w, bz);
  TMmultM_nomem(bz, z, x);
  if (w == NULL) {
    // formula in 'linear' math: X_update = n^-1 * B(Z) * Z
    Mmultscalar_nomem(x, norm);
    return x;
  } else {
    // formula in 'linear' math: X_update = V+ * B(Z) * Z
    TMmultM_nomem(pinv, x, z);
    memcpy(x->elems, z->elems, z->num_elems * sizeof(real_t));
    return x;
  }
}

tm_t *computeV(const tm_t *w) {
  tm_t *v = initTM(w->size);
  unsigned int i, j;
#pragma omp parallel for shared(w, v) private(i) schedule(static)
  for (i = 0; i < w->num_elems; ++i) {
    v->elems[i] = -w->elems[i];
  }
  real_t sum = 0;
#pragma omp parallel for shared(v) private(i, j, sum) schedule(dynamic)
  for (i = 0; i < v->size; ++i) {
    sum = 0;
    for (j = 0; j < v->size; ++j) {
      sum += (j <= i) ? v->elems[(i * (i + 1)) / 2 + j]
                      : v->elems[(j * (j + 1)) / 2 + i];
    }
    v->elems[(i * (i + 1)) / 2 + i] = -sum;
  }
  return v;
}

real_t computesigma(const m_t *x, const tm_t *delta, const tm_t *w) {
  real_t sigma = 0;
  unsigned int i;
  tm_t *dx = calcDeltaMatrix(x, EUCLIDEAN, 2);
  CERROR(dx->size != delta->size, "bad dimension to compute sigma");
  if (w == NULL) {
#pragma omp parallel for shared(delta,dx) private(i) reduction(+:sigma) schedule(static)
    for (i = 0; i < delta->num_elems; ++i) {
      sigma += powf(delta->elems[i] - dx->elems[i], 2);
    }
  } // if we have weights, we have to use a more advanced formula
  else {
#pragma omp parallel for shared(delta,dx,w) private(i) reduction(+:sigma) schedule(static)
    for (i = 0; i < delta->num_elems; ++i) {
      sigma += w->elems[i] * powf(delta->elems[i] - dx->elems[i], 2);
    }
  }
  // cleanup temporary stuff
  freeTM(dx);
  dx = NULL;
  return sigma;
}

real_t computesigma_nomem(const tm_t *dx, const tm_t *delta, const tm_t *w) {
  real_t sigma = 0;
  unsigned int i;
  CERROR(dx->size != delta->size, "bad dimension to compute sigma");
  if (w == NULL) {
#pragma omp parallel for shared(delta,dx) private(i) reduction(+:sigma) schedule(static)
    for (i = 0; i < delta->num_elems; ++i) {
      sigma += powf(delta->elems[i] - dx->elems[i], 2);
    }
  } // if we have weights, we have to use a more advanced formula
  else {
#pragma omp parallel for shared(delta,dx,w) private(i) reduction(+:sigma) schedule(static)
    for (i = 0; i < delta->num_elems; ++i) {
      sigma += w->elems[i] * powf(delta->elems[i] - dx->elems[i], 2);
    }
  }
  return sigma;
}

tm_t *computeBofZ(const tm_t *delta, const tm_t *dz, const tm_t *w) {
  tm_t *b = initTM(delta->size);
  unsigned int i, j;
  if (w == NULL) {
#pragma omp parallel for shared(b, delta, dz, w) private(i, j) schedule(dynamic)
    for (i = 0; i < b->size; ++i) {
      for (j = 0; j < i; ++j) {
        b->elems[(i * (i + 1) / 2) + j] =
            (dz->elems[(i * (i + 1) / 2) + j] != 0)
                ? -(delta->elems[(i * (i + 1) / 2) + j] /
                    dz->elems[(i * (i + 1) / 2) + j])
                : 0.0;
      }
    }
  } else {
#pragma omp parallel for shared(b, delta, dz, w) private(i, j) schedule(dynamic)
    for (i = 0; i < b->size; ++i) {
      for (j = 0; j < i; ++j) {
        b->elems[(i * (i + 1) / 2) + j] =
            (dz->elems[(i * (i + 1) / 2) + j] != 0)
                ? -((w->elems[(i * (i + 1) / 2) + j] *
                     delta->elems[(i * (i + 1) / 2) + j]) /
                    dz->elems[(i * (i + 1) / 2) + j])
                : 0.0;
      }
    }
  }
  // Caution: diagonal entries have to be calculated at last!
  real_t sumcol;
#pragma omp parallel for shared(b) private(i, j, sumcol) schedule(dynamic)
  for (i = 0; i < b->size; ++i) {
    sumcol = 0;
    for (j = 0; j < b->size; ++j) {
      if (j < i)
        sumcol += b->elems[(i * (i + 1) / 2) + j];
      else
        sumcol += b->elems[(j * (j + 1) / 2) + i];
    }
    b->elems[(i * (i + 1) / 2) + i] = -sumcol;
  }
  return b;
}

void computeBofZ_nomem(const tm_t *delta, const tm_t *dz, const tm_t *w,
                       tm_t *b) {
  unsigned int i, j;
  if (w == NULL) {
#pragma omp parallel for shared(b, delta, dz, w) private(i, j) schedule(dynamic)
    for (i = 0; i < b->size; ++i) {
      for (j = 0; j < i; ++j) {
        b->elems[(i * (i + 1) / 2) + j] =
            (dz->elems[(i * (i + 1) / 2) + j] != 0)
                ? -(delta->elems[(i * (i + 1) / 2) + j] /
                    dz->elems[(i * (i + 1) / 2) + j])
                : 0.0;
      }
    }
  } else {
#pragma omp parallel for shared(b, delta, dz, w) private(i, j) schedule(dynamic)
    for (i = 0; i < b->size; ++i) {
      for (j = 0; j < i; ++j) {
        b->elems[(i * (i + 1) / 2) + j] =
            (dz->elems[(i * (i + 1) / 2) + j] != 0)
                ? -((w->elems[(i * (i + 1) / 2) + j] *
                     delta->elems[(i * (i + 1) / 2) + j]) /
                    dz->elems[(i * (i + 1) / 2) + j])
                : 0.0;
      }
    }
  }
  // Caution: diagonal entries have to be calculated at last!
  real_t sumcol;
#pragma omp parallel for shared(b) private(i, j, sumcol) schedule(dynamic)
  for (i = 0; i < b->size; ++i) {
    sumcol = 0;
    for (j = 0; j < b->size; ++j) {
      if (i == j) {

      } else if (j < i) {
        sumcol += b->elems[(i * (i + 1) / 2) + j];
      } else {
        sumcol += b->elems[(j * (j + 1) / 2) + i];
      }
    }
    b->elems[(i * (i + 1) / 2) + i] = -sumcol;
  }
}

m_t *generateRandM(const unsigned int rows, const unsigned int cols,
                   const unsigned int min, const unsigned int max) {
  m_t *mrand = initM(rows, cols);
  // initialize the random generator (never try something parallel with c random
  // functions!)
#if SEEDED == 1
  srand(time(NULL));
#endif
  for (unsigned int i = 0; i < mrand->num_elems; ++i) {
    mrand->elems[i] = min + rand() / (RAND_MAX / (max - min + 1) + 1);
  }
  return mrand;
}

tm_t *convertMtoTM(const m_t *tmp_w) {
  CERROR(tmp_w->rows != tmp_w->cols,
         "weight matrix of type m has wrong dimensions");
  tm_t *w = initTM(tmp_w->rows);
  unsigned int i, j;
#pragma omp parallel for shared(w, tmp_w) private(i, j) schedule(dynamic)
  for (i = 0; i < tmp_w->rows; ++i) {
    for (j = 0; j < i; ++j) {
      w->elems[(i * (i + 1)) / 2 + j] = tmp_w->elems[i * tmp_w->cols + j];
    }
  }
  return w;
}

tm_t *get_test_delta() {
  tm_t *test_delta = initTM(4);
  test_delta->elems[0] = 0;
  test_delta->elems[1] = 5;
  test_delta->elems[2] = 0;
  test_delta->elems[3] = 3;
  test_delta->elems[4] = 2;
  test_delta->elems[5] = 0;
  test_delta->elems[6] = 4;
  test_delta->elems[7] = 2;
  test_delta->elems[8] = 1;
  test_delta->elems[9] = 0;
  return test_delta;
}

m_t *get_test_x() {
  m_t *test_x = initM(4, 2);
  test_x->elems[0] = -0.266;
  test_x->elems[1] = -0.539;
  test_x->elems[2] = 0.451;
  test_x->elems[3] = 0.252;
  test_x->elems[4] = 0.016;
  test_x->elems[5] = -0.238;
  test_x->elems[6] = -0.200;
  test_x->elems[7] = 0.524;
  return test_x;
}
