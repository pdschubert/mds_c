/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/** @file io.c
 *  @brief Implementation der Prototypen aus linalg.h.
 *
 *  In dieser Datei sind die Funktionen zur Realisierung der mathematischen
 *  Operationen auf den Matrixdatenstrukturen implementiert.
 *
 *  @author Philipp D. Schubert
 *  @bug Keine Bugs bekannt.
 */

#include "linalg.h"
#include "m.h"
#include "tm.h"
#include "utils.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

m_t *MmultM(const m_t *A, const m_t *B) {
  CERROR(A->cols != B->rows, "bad dimension for matrix multiplication");

  m_t *C = initM(A->rows, B->cols);
  unsigned int i, j, k;

#pragma omp parallel for shared(A, B, C) private(i, j, k) schedule(static)
  for (i = 0; i < A->rows; ++i) {
    for (k = 0; k < A->cols; ++k) {
      for (j = 0; j < B->cols; ++j) {
        C->elems[i * C->cols + j] +=
            A->elems[i * A->cols + k] * B->elems[k * B->cols + j];
      }
    }
  }
  return C;
}

m_t *MaddM(const m_t *A, const m_t *B) {
  CERROR(A->rows != B->rows || A->cols != B->cols,
         "bad dimension for matrix addition");

  m_t *C = initM(A->rows, A->cols);
  unsigned int i;

#pragma omp parallel for shared(A, B, C) private(i) schedule(static)
  for (i = 0; i < A->rows * A->cols; ++i) {
    C->elems[i] = A->elems[i] + B->elems[i];
  }
  return C;
}

m_t *Mmultscalar(const m_t *A, const real_t x) {
  m_t *C = initM(A->rows, A->cols);
  unsigned int i;

#pragma omp parallel for shared(A, C) private(i) schedule(static)
  for (i = 0; i < A->rows * A->cols; ++i) {
    C->elems[i] = A->elems[i] * x;
  }
  return C;
}

void Mmultscalar_nomem(m_t *A, const real_t x) {
  unsigned int i;
#pragma omp parallel for shared(A) private(i) schedule(static)
  for (i = 0; i < A->rows * A->cols; ++i) {
    A->elems[i] = A->elems[i] * x;
  }
}

m_t *Maddscalar(const m_t *A, const real_t x) {
  m_t *C = initM(A->rows, A->cols);
  unsigned int i;

#pragma omp parallel for shared(A, C) private(i) schedule(static)
  for (i = 0; i < A->rows * A->cols; ++i) {
    C->elems[i] = A->elems[i] + x;
  }
  return C;
}

tm_t *TMaddTM(const tm_t *A, const tm_t *B) {
  CERROR(A->size != B->size, "bad dimension for matrix addition");

  tm_t *C = initTM(A->size);
  unsigned int i;

#pragma omp parallel for shared(A, B, C) private(i) schedule(static)
  for (i = 0; i < A->num_elems; ++i) {
    C->elems[i] = A->elems[i] + B->elems[i];
  }
  return C;
}

tm_t *TMmultscalar(const tm_t *A, const real_t x) {
  tm_t *C = initTM(A->size);
  unsigned int i;

#pragma omp parallel for shared(A, C) private(i) schedule(static)
  for (i = 0; i < A->num_elems; ++i) {
    C->elems[i] = A->elems[i] * x;
  }
  return C;
}

tm_t *TMaddscalar(const tm_t *A, const real_t x) {
  tm_t *C = initTM(A->size);
  unsigned int i;

#pragma omp parallel for shared(A, C) private(i) schedule(static)
  for (i = 0; i < A->num_elems; ++i) {
    C->elems[i] = A->elems[i] + x;
  }
  return C;
}

m_t *MmultTM(const m_t *A, const tm_t *B) {
  CERROR(A->cols != B->size, "bad dimension for matrix multiplication");

  m_t *C = initM(A->rows, B->size);
  unsigned int i, j, k;

#pragma omp parallel for shared(A, B, C) private(i, j, k) schedule(static)
  for (i = 0; i < A->rows; ++i) {
    for (j = 0; j < B->size; ++j) {
      for (k = 0; k < A->cols; ++k) {
        if (k <= j)
          C->elems[i * C->cols + j] +=
              A->elems[i * A->cols + k] * B->elems[k + (j * (j + 1)) / 2];
        else
          C->elems[i * C->cols + j] +=
              A->elems[i * A->cols + k] * B->elems[j + (k * (k + 1)) / 2];
      }
    }
  }
  return C;
}

m_t *TMmultM(const tm_t *A, const m_t *B) {
  CERROR(A->size != B->rows, "bad dimension for matrix multiplication");

  m_t *C = initM(A->size, B->cols);
  unsigned int i, j, k, z;

  real_t abuffer[A->size];
#pragma omp parallel for shared(A, B, C) private(i, j, k, z, abuffer)          \
    schedule(dynamic)
  for (i = 0; i < A->size; ++i) {
    for (z = 0; z < A->size; ++z) {
      abuffer[z] = (z <= i) ? A->elems[z + (i * (i + 1)) / 2]
                            : A->elems[i + (z * (z + 1)) / 2];
    }
    for (j = 0; j < B->cols; ++j) {
      for (k = 0; k < A->size; ++k) {
        C->elems[i * C->cols + j] += abuffer[k] * B->elems[k * B->cols + j];
      }
    }
  }
  return C;
}

void TMmultM_nomem(const tm_t *A, const m_t *B, m_t *C) {
  CERROR(A->size != B->rows, "bad dimension for matrix multiplication");
  memset(C->elems, 0, C->num_elems * sizeof(real_t));
  unsigned int i, j, k, z;
  real_t abuffer[A->size];
#pragma omp parallel for shared(A, B, C) private(i, j, k, z, abuffer)          \
    schedule(dynamic)
  for (i = 0; i < A->size; ++i) {
    for (z = 0; z < A->size; ++z) {
      abuffer[z] = (z <= i) ? A->elems[z + (i * (i + 1)) / 2]
                            : A->elems[i + (z * (z + 1)) / 2];
    }
    for (j = 0; j < B->cols; ++j) {
      for (k = 0; k < A->size; ++k) {
        C->elems[i * C->cols + j] += abuffer[k] * B->elems[k * B->cols + j];
      }
    }
  }
}
