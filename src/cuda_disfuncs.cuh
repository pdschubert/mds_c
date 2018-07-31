/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file cuda_disfuncs.cuh
 * @brief Funktionsprototypen für die Berechnung der Distanzmatrix d.
 *
 * Diese Datei enthält die Funktionsprototypen für die Berechnung von
 * Distanzmatrizen, sowie eine Testfunktion, die eine korrekte Funktionsweise
 * eines Kernel-Aufrufes feststellt.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef CUDA_DISFUNCS_CUH_
#define CUDA_DISFUNCS_CUH_

#include "m.h"
#include "tm.h"
#include "utils.h"

/**
 * @brief Führt einen Testkernel aus.
 *
 * Es wird ein Testkernel ausgeführt, der eine Ausgabe auf der Kommandozeile
 * liefert. Bei einer Fehlermeldung ist etwas auf dem System falsch
 * konfiguriert.
 */
extern void cuda_disfuncs_test();

/**
 * @brief Berechnet die Distanzmatrix nach der euklidischen Metrik.
 *
 * Es wird die euklidsche Distanzmatrix einer gegebenen Datenmatrix berechnet.
 * Die Einträge der Matrix stellen die paarweisen Distanzen dar.
 *
 * @param data eine beliebige Datenmatrix
 * @return Dreiecksmatrix, die die paarweisen Distanzen enthält
 */
extern tm_t *cuda_distanceMatrix(const m_t *data);

/**
 * @brief Berechnet die Distanzmatrix nach der euklidischen Metrik ohne eigenes
 * Speichermanagement.
 *
 * Die Funktion berechnet die euklidische Distanzmatrix einer gegebenen
 * Datematrix. Es wird von der Funktion dabei kein eigener Speicher alloziert.
 *
 * @warning Es dürfen ausschließlich Device Pointer übergeben werden!
 * @param dev_x Device Pointer der Datenmatrix
 * @param xrows Anzahl der Zeiler der Datenmatrix
 * @param xcols Anzahl der Spalten der Datenmatrix
 * @param dev_d Device Pointer an dessen Stelle die Distanzmatrix geschrieben
 * werden soll
 *
 */
extern void cuda_distanceMatrix_nomem(const real_t *dev_x, unsigned int xrows,
                                      unsigned int xcols, real_t *dev_d);

#endif /* CUDA_DISFUNCS_CUH_ */
