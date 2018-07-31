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
 * @brief Enthält die Datenstruktur zur Speicherung allgemeiner Matrizen, sowie
 * grundlegende Funktionen.
 *
 * Diese Datei enthält die Definitionen für die Datenstrukturen zur Speicherung
 * allgemeiner Matrizen, sowie einige elementare Funktionen, um mit den
 * Datenstrukturen umzugehen.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef M_H_
#define M_H_

#include "utils.h"

/**
 * @brief Datenstruktur zur Speicherung von Matrizen.
 */
typedef struct mat {
  unsigned int rows;      /**< Anzahl der Zeilen */
  unsigned int cols;      /**< Anzahl der Spalten */
  unsigned int num_elems; /**< Anzahl der Matrixelemente */
  real_t
      *elems; /**< Adresse des Datenblocks, an dem die Matrixelemente liegen */
} m_t;

/**
 * @brief Initialisiert eine allgemeine Matrix.
 *
 * Es wird eine allgemeine Matrix vorgegebener Größe initalisiert, Speicher für
 * diese angefordert und der Adresse der initialisierten Matrix zurückgeliefert.
 *
 * @param rows Anzahl der Zeilen
 * @param cols Anzahl der Zeilen
 * @return Adresse der initialisierten Matrix
 */
m_t *initM(const unsigned int rows, const unsigned int cols);

/**
 * @brief Gibt den Speicher einer allgemeinen Matrix frei.
 *
 * Der Speicher der übergebenen allgemeinen Matrix wird freigegeben.
 *
 * @param m Matrix zum Freigeben
 */
void freeM(m_t *m);

/**
 * @brief Gibt eine Matrix auf der Kommandozeile aus.
 *
 * Die Funktion gibt die übergebene Matrix formatiert auf der Kommandozeile
 * aus.
 *
 * @param m Matrix die ausgegeben werden soll
 */
void printM(const m_t *m);

#endif /* M_H_ */
