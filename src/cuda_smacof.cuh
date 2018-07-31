/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file cuda_smacof.cuh
 * @brief Funktionsprototyp für die Funktion, die den SMACOF Algorithmus in
 * MEMPOOL Variante durchführt.
 *
 * Diese Datei enthält die Funktionsprototypen für die Funktion zur Berechnung
 * des SMACOF Algorithmus. Dies ist die CUDA Variante, die einen Memory Pool
 * benutzt. D.h. es werden beim Start des Algorithmus alle zur Berechnung
 * benötigten Matrix alloziert und mit diesen wird anschließend gearbeitet. Zur
 * Verarbeitung müssen dann die SMACOF Hilfsfunktionen in der 'nomem' bzw.
 * 'noret' Variante verwendet werden.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef CUDA_SMACOF_CUH_
#define CUDA_SMACOF_CUH_

#include "m.h"
#include "tm.h"
#include "utils.h"

/**
 * @brief Berechnet den SMACOF Algorithmus und gibt die Lösungskonfiguration
 * zurück.
 *
 * Es wird der SMACOF Algorithmus berechnet, wobei ein Memory Pool verwendet
 * wird.
 *
 * @param delta Unähnlichkeitsmatrix
 * @param w Gewichtsmatrix
 * @param pinv Pseudoinverses V+
 * @param epsilon Schwellwert, für vorzeitigen Abbruch
 * @param dim Ergebnisdimension
 * @param maxiter Maximale Zahl an Iterationen, die benutzt werden
 * @param print Programm gibt sich gesprächig bzw. nicht
 * @return Lösungskonfiguration
 */
extern m_t *cuda_computeSMACOF_mempool(const tm_t *delta, const tm_t *w,
                                       const tm_t *pinv, const real_t epsilon,
                                       const unsigned int dim,
                                       const unsigned int maxiter,
                                       const unsigned int print);

#endif /* CUDA_SMACOF_CUH_ */
