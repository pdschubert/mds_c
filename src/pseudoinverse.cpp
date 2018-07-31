/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file pseudoinverse.cpp
 * @brief Enthält die Implementation der Datei pseudoinverse.hh.
 *
 * Hier wird der Prototyp zur Berechnung eines Pseudoinversen einer
 * Dreiecksmatrix implementiert.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#include <armadillo>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>

extern "C" {
#include "m.h"
#include "tm.h"
#include "utils.h"
}

extern "C" tm_t *pseudoinverse(const tm_t *t) {
  m_t *m = initM(t->size, t->size);
  unsigned int i, j;
#pragma omp parallel for shared(t, m) private(i, j) schedule(dynamic)
  for (i = 0; i < t->size; ++i) {
    for (j = 0; j < t->size; ++j) {
      m->elems[i * t->size + j] = (j <= i) ? t->elems[(i * (i + 1)) / 2 + j]
                                           : t->elems[(j * (j + 1)) / 2 + i];
    }
  }
  // check if we should use single or double precision and use adequate
  // armadillo data structures
  if (typeid(real_t).name() == typeid(float).name()) {
    arma::fmat f_am =
        arma::fmat((float *)m->elems, t->size, t->size, false, false);
    arma::fmat f_inverse = arma::pinv(f_am);
    // caution: override data memory of m
    memcpy(m->elems, f_inverse.memptr(), t->size * t->size * sizeof(real_t));
  } else {
    arma::mat d_am =
        arma::mat((double *)m->elems, t->size, t->size, false, false);
    arma::mat d_inverse = arma::pinv(d_am);
    // caution: override data memory of m
    memcpy(m->elems, d_inverse.memptr(), t->size * t->size * sizeof(real_t));
  }
  tm_t *result = initTM(t->size);
#pragma omp parallel for shared(result, m) private(i, j) schedule(dynamic)
  for (i = 0; i < t->size; ++i) {
    for (j = 0; j <= i; ++j) {
      result->elems[(i * (i + 1)) / 2 + j] = m->elems[i * t->size + j];
    }
  }
  free(m);
  return result;
}
