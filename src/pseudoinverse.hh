/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file pseudoinverse.hh
 * @brief Enthält den Prototypen zur Berechnung eines Pseudoinversen einer
 * Dreiecksmatrix.
 *
 * Hier wird der Prototyp zur Berechnung eines Pseudoinversen einer
 * Dreiecksmatrix deklariert. Diese wird benötigt, um das Pseudoinverse V+ zur
 * Matrix V zu berechnen.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef PSEUDOINVERSE_H_
#define PSEUDOINVERSE_H_

#include "tm.h"

/**
 * @brief Berechnet das Pseudoinverse der übergebenen Dreiecksmatrix.
 *
 * Die Funktion berechnet das Pseudoinverse einer Dreiecksmatrix mit Hilfe der
 * Armadillo C++ Bibliothek. Es wurde dabei auf eine besonders effiziente
 * Typumwandlung von den hier eingesetzten Datenstrukturen in die
 * Armadillo-Eigenen Datenstrukturen wertgelegt. Die durchgeführten
 * Typumwandlungen besitzen nahezu keinen Mehrwaufwand, da ausschließlich
 * Pointer übergeben werden, um den Datenteil der Matrizen hin und her
 * zureichen. Das Pseudoinverse einer quadratischen, symmetrischen Matrix ist
 * wiederum quadratisch und symmetrisch, weswegen eine Dreiecksmatrix
 * zurückgegeben werden kann. Innerhalb der Funktion muss die aktuell verwendete
 * Darstellung von Gleitkommazahlen ermittelt werden, um in Armadillo das
 * entsprechende Pendant zu wählen.
 *
 * @warning Die Funktion ist für Kompatibilität mit Armadillo in C++
 * implementiert.
 * @param t Dreiecksmatrix
 * @return Pseudoinverse von t als Dreiecksmatrix
 */
extern tm_t *pseudoinverse(const tm_t *t);

#endif /* PSEUDOINVERSE_H_ */
