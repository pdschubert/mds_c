/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file tm.c
 * @brief Enthält die Implementation der Prototypen aus tm.h.
 *
 * Diese Datei enthält die Implementation der Prototypen der elementaren
 * Funktionen zur Verwaltung von Dreiecksmatrizen.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#include "tm.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

tm_t *initTM(const unsigned int size) {
  tm_t *m;
  if ((m = (tm_t *)malloc(sizeof(tm_t) + ((size * (size + 1)) / 2) *
                                             sizeof(real_t))) == NULL) {
    perror("memory allocation failed");
    HEREANDNOW;
    exit(-1);
  }
  m->size = size;
  m->num_elems = (size * (size + 1)) / 2;
  m->elems = (real_t *)(&m->elems + 1);
  memset(m->elems, 0, m->num_elems * sizeof(real_t));
  return m;
}

void freeTM(tm_t *m) { free(m); }

void printTM(const tm_t *m) {
  if (m == NULL) {
    puts("empty matrix of type tm");
    return;
  }
  for (unsigned int i = 0; i < m->size; i++) {
    for (unsigned int j = 0; j <= i; j++) {
      printf("%f ", m->elems[j + (i * (i + 1)) / 2]);
    }
    puts("");
  }
}
