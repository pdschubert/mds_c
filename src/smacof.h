/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file smacof.h
 * @brief Enthält die Prototypen der Funktionen zur Berechnung des SMACOF
 * Algorithmus.
 *
 * Die Datei enthält die Prototypen zur Berechnung des SMACOF Algorithmus in C
 * in der Standard- und in der Memory-Pool-Version.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef SMACOF_H_
#define SMACOF_H_

#include "m.h"
#include "tm.h"
#include "utils.h"

/**
 * @brief Berechnet die Lösung des MDS Problems mittels SMACOF Verfahren.
 *
 * Es wird eine Lösung des Problems der multidimensionalen Skalierung mit Hilfe
 * des iterativen SMACOF Verfahrens berechnet. Die Berechnungen werden mit einem
 * Standarddesign aller Funktionen durchgeführt.
 *
 * @param delta Unähnlichkeitsmatrix
 * @param w Gewichtsmatrix
 * @param pinv Pseudoinverse V+
 * @param epsilon Schwellwert für vorzeitigen Abbruch
 * @param dim Dimensionalität der Lösung
 * @param maxiter Maximale Anzahl an Iterationen die benutzt werden soll
 * @param print Gesprächigkeit des Algorithmus
 * @result Lösungskonfiguration
 */
m_t *computeSMACOF(const tm_t *delta, const tm_t *w, const tm_t *pinv,
                   const real_t epsilon, const unsigned int dim,
                   const unsigned int maxiter, const unsigned int print);

/**
 * @brief Berechnet die Lösung des MDS Problems mittels SMACOF Verfahren.
 *
 * Es wird eine Lösung des Problems der multidimensionalen Skalierung mit Hilfe
 * des iterativen SMACOF Verfahrens berechnet. Die Berechnungen werden mit einem
 * Memory-Pool und einem speziellen Design einiger Funktionen durchgeführt.
 * Zu Beginn des Algorithmus werden alle benötigten Speicherstrukturen alloziert
 * und es wird im weiteren Algorithmus ausschließlich mit diesen gearbeitet.
 * Dies bringt einen großen Geschwindgkeitsvorteil, allerdings auf Kosten der
 * Lesbarkeit des Programms.
 *
 * @param delta Unähnlichkeitsmatrix
 * @param w Gewichtsmatrix
 * @param pinv Pseudoinverse V+
 * @param epsilon Schwellwert für vorzeitigen Abbruch
 * @param dim Dimensionalität der Lösung
 * @param maxiter Maximale Anzahl an Iterationen die benutzt werden soll
 * @param print Gesprächigkeit des Algorithmus
 * @result Lösungskonfiguration
 */
m_t *computeSMACOF_mempool(const tm_t *delta, const tm_t *w, const tm_t *pinv,
                           const real_t epsilon, const unsigned int dim,
                           const unsigned int maxiter,
                           const unsigned int print);

/**
 * @brief Funktion zum Testen des CUDA Funktionssatzes.
 */
void test_cuda_functions();

#endif /* SMACOF_H_ */
