/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/** @file disfuncs.c
 *  @brief Implementation der Prototypen aus disfuncs.h.
 *
 *  In dieser Datei sind die Funktionen zur Berechnung der Unähnlichkeitsmatrix
 * Delta, sowie die zur Berechnung einer (euklidischen) Distanzmatrix (ohne
 * eigenes Speichermanagment).
 *
 *  @author Philipp D. Schubert
 *  @bug Keine Bugs bekannt.
 */

#include "disfuncs.h"
#include "m.h"
#include "tm.h"
#include "utils.h"
#include <omp.h>
#include <stdlib.h>
#include <string.h>

tm_t *calcDeltaMatrix(const m_t *data, const enum delta d, const int p) {
  tm_t *delta = initTM(data->rows);
  unsigned int i, j;
  switch (d) {
  case EUCLIDEAN:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            euclidean(&data->elems[i * data->cols],
                      &data->elems[j * data->cols], data->cols);
      }
    }
    break;
  case CITYBLOCK:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            cityblock(&data->elems[i * data->cols],
                      &data->elems[j * data->cols], data->cols);
      }
    }
    break;
  case MINKOWSKI:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            minkowski(&data->elems[i * data->cols],
                      &data->elems[j * data->cols], p, data->cols);
      }
    }
    break;
  case CANBERRA:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            canberra(&data->elems[i * data->cols], &data->elems[j * data->cols],
                     data->cols);
      }
    }
    break;
  case DIVERGENCE:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            divergence(&data->elems[i * data->cols],
                       &data->elems[j * data->cols], data->cols);
      }
    }
    break;
  case BRAYCURTIS:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            braycurtis(&data->elems[i * data->cols],
                       &data->elems[j * data->cols], data->cols);
      }
    }
    break;
  case SOERGEL:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            soergel(&data->elems[i * data->cols], &data->elems[j * data->cols],
                    data->cols);
      }
    }
    break;
  case BAHATTACHARYYA:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            bahattacharyya(&data->elems[i * data->cols],
                           &data->elems[j * data->cols], data->cols);
      }
    }
    break;
  case WAVEHEDGES:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            wavehedges(&data->elems[i * data->cols],
                       &data->elems[j * data->cols], data->cols);
      }
    }
    break;
  case ANGULARSEPERATION:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            angularseperation(&data->elems[i * data->cols],
                              &data->elems[j * data->cols], data->cols);
      }
    }
    break;
  case CORRELATION:
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j < i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] =
            correlation(&data->elems[i * data->cols],
                        &data->elems[j * data->cols], data->cols);
      }
    }
    break;
  case NONE:
    CERROR(data->rows != data->cols, "input is no valid dissimilarity matrix");
#pragma omp parallel for shared(delta) private(i, j) schedule(dynamic)
    for (i = 0; i < data->rows; ++i) {
      for (j = 0; j <= i; ++j) {
        delta->elems[(i * (i + 1)) / 2 + j] = data->elems[i * data->cols + j];
      }
    }
    break;
  default:
    perror("unsupported measure of DELTA");
    HEREANDNOW;
    exit(-1);
    break;
  }
  return delta;
}

void calcDistanceMatrix_nomem(const m_t *x, tm_t *dx) {
  unsigned int i, j;
#pragma omp parallel for shared(x, dx) private(i, j) schedule(dynamic)
  for (i = 0; i < x->rows; ++i) {
    for (j = 0; j < i; ++j) {
      dx->elems[(i * (i + 1)) / 2 + j] =
          euclidean(&x->elems[i * x->cols], &x->elems[j * x->cols], x->cols);
    }
  }
}
