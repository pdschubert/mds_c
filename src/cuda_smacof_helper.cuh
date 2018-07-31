/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file cuda_smacof_helper.cuh
 * @brief Funktionsprototypen der für den SMACOF Algorithmus benötigten
 * Funktionen.
 *
 * Diese Datei enthält die Funktionsprototypen für die benötigten Funktionen,
 * zur Formulierung des SMACOF Algorithmus.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef CUDA_SMACOF_HELPER_CUH_
#define CUDA_SMACOF_HELPER_CUH_

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
extern void cuda_smacof_helper_test();

/**
 * @brief Führt eine Guttman-Transformation auf die gewohnte Art und Weise
 * durch.
 *
 * Es wird eine Guttman-Transformation durchgeführt. Das errechnete, update
 * Konfiguration wird anschließend zurückgegeben. Wird w = NULL übergeben,
 * so kann die einfachere Updateformel verwendet werden, sind Gewichte
 * angegeben, d.h. w != NULL, so wird die erweiterte Updateformel verwendet.
 *
 * @param delta Unähnlichkeitsmatrix
 * @param z aktuelle Konfiguration
 * @param w Gewichtsmatrix
 * @param pinv Pseudoinverses V+
 * @return Konfigurationsupdate
 */
extern m_t *cuda_guttmanTransformation(const tm_t *delta, const m_t *z,
                                       const tm_t *w, const tm_t *pinv);

/**
 * @brief Führt eine Guttman-Transformation ohne eigenes Speichermanagement
 * durch.
 *
 * Es wird eine Guttman-Transformation durchgeführt, ohne das diese Funktion
 * ein Speichermanagement durchführt. Wird w = NULL übergeben,
 * so kann die einfachere Updateformel verwendet werden, sind Gewichte
 * angegeben, d.h. w != NULL, so wird die erweiterte Updateformel verwendet.
 *
 * @warning Es dürfen ausschließlich Device Pointer übergeben werden
 * @param dev_delta Array der Unähnlichkeiten Delta
 * @param delta_size Größe der Delta Matrix
 * @param delta_num_elems Länge des Arrays dev_delta
 * @param dev_z Array mit der aktuellen Konfiguration
 * @param zrows Anzahl der Zeilen von z
 * @param zcols Anzahl der Spalten von z
 * @param dev_x Array an welches das Ergebnis geschrieben wird
 * @param dev_w Array mit der Gewichtung
 * @param dev_pinv Array mit Pseudoinversem V+
 * @param dev_d Array mit Distanzmatrix
 * @param dev_bz Array mit B Matrix
 */
extern real_t *cuda_guttmanTransformation_nomem(
    const real_t *dev_delta, unsigned int delta_size,
    unsigned int delta_num_elems, real_t *dev_z, unsigned int zrows,
    unsigned int zcols, real_t *dev_x, const real_t *dev_w,
    const real_t *dev_pinv, const real_t *dev_d, real_t *dev_bz);

/**
 * @brief Berechnet aus der Gewichtsmatrix die Matrix V.
 *
 * Es wird aus der übergebenen Gewichtsmatrix W, die Matrix V
 * berechnet. Diese enthält die negierten Einträge von W, sowie
 * die negierten Spaltensummen auf der Hauptdiagonalen.
 *
 * @param w Gewichtsmatrix
 * @return Matrix V
 */
extern tm_t *cuda_computeV(const tm_t *w);

/**
 * @brief Berechnet den Stresswert sigma.
 *
 * Berechnet den Stresswert sigma aus der übergebenen Konfiguration
 * unter Berücksichtigung der Gewichtsmatrix.
 *
 * @param x Konfiguration
 * @param delta Unähnlichkeitsmatrix
 * @param w Gewichtsmatrix
 * @return sigma Stress der Konfiguration
 */
extern real_t cuda_computesigma(const m_t *x, const tm_t *delta, const tm_t *w);

/**
 * @brief Berechnet den Stresswert sigma ohne eigenes Speichermanagement.
 *
 * Es wird der Stresswert einer Konfiguration berechnet ohne das eine
 * Speichermanagement durch die Funktion durchgeführt wird.
 *
 * @warning Es dürfen ausschließlich Device Pointer übergeben werden
 * @param dev_d Array der Distanzmatrix
 * @param dnum_elems Länge des Arrays dev_d
 * @param dev_delta Array der Unähnlichkeitsmatrix
 * @param dev_w Array der Gewichtsmatrix
 * @return sigma Stress der Konfiguration
 */
extern real_t cuda_computesigma_nomem(const real_t *dev_d,
                                      const unsigned int dnum_elems,
                                      const real_t *dev_delta,
                                      const real_t *dev_w);

/**
 * @brief Berechnet die B Matrix aus der aktuellen Distanzmatrix.
 *
 * Es wird die Matrix B aus der aktuellen Distanzmatrix und der
 * Unähnlichkeitsmatrix berechnet. Ist w != NULL, also eine Gewichtsmatrix
 * angegeben, so wird diese bei der Berechnung von B berücksichtigt. B enthält
 * die negativen (gewichteten) Quotienten der Einträge der Unähnlichkeits- und
 * Distanzmatrix oder 0.0, falls durch Nenner 0.0 ist.
 *
 * @param delta Unähnlichkeitsmatrix
 * @param dz Distanzmatrix der aktuellen Konfiguration
 * @param w Gewichtsmatrix
 * @return Matrix B
 */
extern tm_t *cuda_computeBofZ(const tm_t *delta, const tm_t *dz, const tm_t *w);

/**
 * @brief Berechnet die B Matrix aus der aktuellen Distanzmatrix ohne
 * Speichermanagement.
 *
 * Es wird die B Matrix aus der aktuellen Distanzmatrix und der
 * Unähnlichkeitsmatrix berechnet. Evtl. (dev_w != NULL) erfolgt eine
 * Gewichtung. Die Funktion führt kein eigenes Speichermanagement durch.
 *
 * @warning Es dürfen nur Device Pointer übergeben werden
 * @param dev_delta Array der Unähnlichkeitsmatrix
 * @param delta_size Größe der Unähnlichkeitsmatrix
 * @param delta_num_elems Länge des Arrays dev_delta
 * @param dev_d Distanzmatrix der aktuellen Konfiguration
 * @param dev_w Array der Gewichtsmatrix
 * @param dev_bz Stelle an die das Ergebnis geschrieben werden soll
 */
extern void cuda_computeBofZ_nomem(const real_t *dev_delta,
                                   const unsigned int delta_size,
                                   const unsigned int delta_num_elems,
                                   const real_t *dev_d, const real_t *dev_w,
                                   real_t *dev_bz);

#endif /* CUDA_SMACOF_HELPER_CUH_ */
