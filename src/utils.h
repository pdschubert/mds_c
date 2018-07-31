/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file utils.h
 * @brief Definiert einige nützliche Makros, sowie den zu benutzenden Typ zur
 * Gleitkommadarstellung.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>

/**
 * Gibt die aktuelle Datei, Zeile und Funktion aus, in der das Makro erreicht
 * wird.
 */
#define HEREANDNOW                                                             \
  fprintf(stderr, "error in :: file: %s line: %d function: %s\n", __FILE__,    \
          __LINE__, __func__);

/**
 * Gibt bei positiver Auswertung eine entsprechende Fehlermeldung aus.
 */
#define CERROR(BOOL, MESSAGE)                                                  \
  if (BOOL) {                                                                  \
    perror(MESSAGE);                                                           \
    HEREANDNOW;                                                                \
    exit(-1);                                                                  \
  }

/**
 * Datentyp der zur Repräsentation von Gleitkommazahlen benutzt werden soll.
 */
typedef float real_t;

/**
 * Die Aufzählung der möglichen Unähnlichkeitsfunktionen.
 */
enum delta {
  EUCLIDEAN = 0,
  CITYBLOCK = 1,
  MINKOWSKI = 2,
  CORRELATION = 3,
  ANGULARSEPERATION = 4,
  WAVEHEDGES = 5,
  BAHATTACHARYYA = 6,
  SOERGEL = 7,
  BRAYCURTIS = 8,
  DIVERGENCE = 9,
  CANBERRA = 10,
  NONE = 11
};

#endif /* UTILS_H_ */
