/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file m.h
 * @brief Implementiert die Prototypen aus m.h.
 *
 * Diese Datei enthält die Implementationen der Prototypen aus m.h zur
 * Verwaltung von Matrizen.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#include "m.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

m_t *initM(const unsigned int rows, const unsigned int cols) {
  m_t *m;
  if ((m = (m_t *)malloc(sizeof(m_t) + rows * cols * sizeof(real_t))) == NULL) {
    perror("memory allocation failed");
    HEREANDNOW;
    exit(-1);
  }
  m->rows = rows;
  m->cols = cols;
  m->num_elems = rows * cols;
  m->elems = (real_t *)(&m->elems + 1);
  memset(m->elems, 0, m->num_elems * sizeof(real_t));
  return m;
}

void freeM(m_t *m) { free(m); }

void printM(const m_t *m) {
  if (m == NULL) {
    puts("empty matrix of type m");
    return;
  }
  for (unsigned int i = 0; i < m->rows; i++) {
    for (unsigned int j = 0; j < m->cols; j++) {
      printf("%f ", m->elems[i * m->cols + j]);
    }
    puts("");
  }
}
