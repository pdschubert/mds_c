/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file tm.h
 * @brief Enthält die Datenstruktur zur Speicherung von Dreiecksmatrizen, sowie
 * grundlegende Funktionen.
 *
 * Diese Datei enthält die Definitionen für die Datenstrukturen zur Speicherung
 * von Dreiecksmatrizen, sowie elementare Funktionen, um diese zu verwalten.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef TM_H_
#define TM_H_

#include "utils.h"

/**
 * @brief Datenstruktur zur Speicherung von quadratischen, symmetrischen
 * Matrizen.
 */
typedef struct trimat {
  unsigned int size; /**< Größe der Matrix, entspricht Anzahl der Zeilen =
                        Anzahl der Spalten */
  unsigned int num_elems; /**< Anzahl der Matrixeinträge */
  real_t
      *elems; /**< Adresse des Datenblocks, an dem die Matrixelemente liegen */
} tm_t;

/**
 * @brief Initialisiert eine Dreiecksmatrix.
 *
 * Es wird eine Matrix gegebener Größe initialisiert, Speicher angefordert und
 * die Adresse der Matrix zurückgegeben.
 *
 * @param size Anzahl der Zeilen = Anzahl der Spalten
 */
tm_t *initTM(const unsigned int size);

/**
 * @brief Gibt den Speicher einer Dreiecksmatrix frei.
 *
 * Der Speicher der übergebenen Dreiecksmatrix wird freigegeben.
 *
 * @param m Matrix zum Freigeben
 */
void freeTM(tm_t *m);

/**
 * @brief Gibt eine Matrix auf der Kommandozeile aus.
 *
 * Die Funktion gibt die übergebene Matrix formatiert auf der Kommandozeile
 * aus.
 *
 * @param m Matrix die ausgegeben werden soll
 */
void printTM(const tm_t *m);

#endif /* TM_H_ */
