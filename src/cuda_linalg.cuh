/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file cuda_linalg.cuh
 * @brief Funktionsprototypen für die mathematischen Operationen auf den
 * Matrizen.
 *
 * Diese Datei enthält die Funktionsprototypen für die üblichen mathematischen
 * Operationen auf den allgemeinen, sowie den Dreiecksmatrizen. Darüberhinaus
 * auch für einen Kernel-Test.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef CUDA_LINALG_CUH_
#define CUDA_LINALG_CUH_

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
extern void cuda_linalg_test();

/**
 * @brief Wrapper-Funktion für eine Multiplikation zwischen zwei allgemeinen
 * Matrizen.
 *
 * Die Wrapper-Funktion, die eine allgemeine Matrixmultiplikation durchführt.
 * In dieser werden die benötigten Vorbereitungen getroffen und anschließend ein
 * CUDA Kernel zur Matrixmultiplikation aufgerufen, dir die Multiplikation
 * realisiert.
 *
 * @param a eine allgemeine Matrix
 * @param b eine allgemeine Matrix
 * @return das Matrixprodukt von a und b
 *
 */
extern m_t *cuda_MmultM(const m_t *a, const m_t *b);

/**
 * @brief Addiert zwei allgemeine Matrizen gleicher Dimension.
 *
 * Es werden zwei allgemeine Matrizen gleicher Dimensionen addiert. Dazu werden
 * die Vorbereitungen getroffen und die Addition mittels eines CUDA Kernels
 * durchgeführt.
 *
 * @param a allgemeine Matrix
 * @param b allgemeine Matrix
 * @return Matrixsumme von a und b
 */
extern m_t *cuda_MaddM(const m_t *a, const m_t *b);

/**
 * @brief Berechnet das Vielfache einer allgemeinen Matrix.
 *
 * Es wird das Vielfache einer allgemeinen Matrix berechnet. Die Berechnungen
 * werden mittels eines CUDA Kernels realisiert.
 *
 * @param a allgemeine Matrix
 * @param x Vielfaches, das berechnet werden soll
 * @return Matrixvielfaches
 */
extern m_t *cuda_Mmultscalar(const m_t *a, const real_t x);

/**
 * @brief Berechnet das Vielfache einer Matrix ohne Speichermanagement.
 *
 * Es wird das Vielfache einer Matrix berechnet ohne das Speicher angefordert
 * wird. Das Ergebnis wird direkt in die übergebenen Matrix zurückgeschrieben.
 *
 * @warning es dürfen ausschließlich Device Pointer übergeben werden
 * @param m eindimensionales Array, dass die Matrix enthält
 * @param x Vielfaches das berechnet werden soll
 * @param len Länge des Arrays m
 */
extern void cuda_Mmultscalar_nomem(real_t *m, const real_t x,
                                   const unsigned int len);

/**
 * @brief Addiert einen skalaren Wert auf jeden Matrixeintrag.
 *
 * Es wird ein skalarer Wert auf jeden Matrixeintrag addiert.
 *
 * @param m allgemeine Matrix
 * @param x Wert der auf die Matrix addiert werden soll
 * @return Ergebnismatrix
 */
extern m_t *cuda_Maddscalar(const m_t *m, const real_t x);

/**
 * @brief Addiert zwei Dreiecksmatrizen gleicher Dimension.
 *
 * Es werden zwei Dreiecksmatrizen der gleichen Dimension addiert.
 * Die Addition wird mittels CUDA Kernel realisiert.
 *
 * @param a Dreiecksmatrix
 * @param b Dreiecksmatrix
 * @return Ergebnismatrix
 */
extern tm_t *cuda_TMaddTM(const tm_t *a, const tm_t *b);

/**
 * @brief Berechnet das Vielfache einer Dreiecksmatrix.
 *
 * Es wird das Vielfache einer Dreiecksmatrix berechnet. Die Multiplikationen
 * der Eintrage werden auf der GPU durchgeführt.
 *
 * @param m Dreiecksmatrix
 * @param x Vielfache
 * @return Matrixvielfaches
 */
extern tm_t *cuda_TMmultscalar(const tm_t *m, const real_t x);

/**
 * @brief Addiert einen skalaren Wert auf alle Matrixeinträge.
 *
 * Es wird ein Skalar auf alle Einträge der Dreiecksmatrix addiert.
 *
 * @param m Dreiecksmatrix
 * @param x Skalar
 * @return Ergebnismatrix
 */
extern tm_t *cuda_TMaddscalar(const tm_t *m, const real_t x);

/**
 * @brief Berechnet das Matrixprodukt aus einer Dreiecks- und einer allgemeinen
 * Matrix.
 *
 * Es wird das Matrixprodukt aus einer Dreiecksmatrix und einer allgemeinen
 * Matrix berechnet. Die Berechnung wird mit einem GPU Kernel realisiert.
 * Die Funktion liefert eine allgemeine Matrix zurück.
 *
 * @param a Dreiecksmatrix
 * @param b allgemeine Matrix
 * @return allgemeine Matrix
 */
extern m_t *cuda_TMmultM(const tm_t *a, const m_t *b);

/**
 * @brief Berechnet das Matrixprodukt aus einer Dreiecks- und einer allgemeinen
 * Matrix, ohne das Speichermanagement dafür zu übernehmen.
 *
 * Es wird das Matrixprodukt einer Dreiecks- und allgemeinen Matrix gebildet,
 * wobei das Speichermanagement nicht von dieser Funktion übernommen wird. Der
 * Funktion dürfen dabei ausschließlich Device Pointer übergeben werden.
 *
 * @warning Es dürfen nur Device Pointer übergeben werden
 * @param dev_a Array der Dreiecksmatrix
 * @param asize Größe der Dreiecksmatrix
 * @param dev_b Array der allgemeinen Matrix
 * @param brows Zeilen der Matrix b
 * @param bcols Spalten der Matrix b
 * @param dev_c Array in welches die Ergebnismatrix geschrieben werden soll
 */
extern void cuda_TMmultM_nomem(const real_t *dev_a, const unsigned int asize,
                               const real_t *dev_b, const unsigned int brows,
                               const unsigned int bcols, real_t *dev_c);

#endif /* CUDA_LINALG_CUH_ */
