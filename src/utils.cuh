/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file utils.cuh
 * @brief Definiert ein nützliches Makro, sowie zugehörige Funktion, zur
 * Überprüfung von CUDA API Funktionen.
 *
 * Hier wird ein nützliches Makro, inklusive der benötigten zugehörigen Funktion
 * definiert, um die Fehlercodes von Funktionen der CUDA API zu überprüfen und
 * im Fehlerfall eine entsprechende Meldung auf der Kommandozeile auszugeben.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef UTILS_CUH_
#define UTILS_CUH_

#include "cuda.h"
#include "cuda_runtime_api.h"
#include "driver_types.h"
#include "stdio.h"

/**
 * @brief Prüft auf Fehlercodes und gibt im Fehlerfall eine entsprechende
 * Meldung aus.
 *
 * Die Funktion überprüft den durch das Makro erhaltenen Rückgabewert einer CUDA
 * API Funkion. Falls einer Fehlerwert festgestellt wird, so wird eine
 * Fehlermeldung mit Informationen zum Ort des Fehler auf der Kommandozeile
 * ausgegeben.
 *
 * @param err Fehlercode
 * @param file aktuelle Datei
 * @param line aktuelle Zeile
 */
static void HandleError(cudaError_t err, const char *file, int line) {
  if (err != cudaSuccess) {
    printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
    exit(EXIT_FAILURE);
  }
}

/**
 * Makro mit dem CUDA API Funktionen ummantelt werden können, um die
 * Rückgabewerte auf Fehlercodes zu überprüfen.
 */
#define CUERROR(err) (HandleError(err, __FILE__, __LINE__))

#endif /* UTILS_CUH_ */
