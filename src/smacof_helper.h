/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file smacof_helper.h
 * @brief Funktionsprototypen der für den SMACOF Algorithmus benötigten
 * Funktionen.
 *
 * Diese Datei enthält die Funktionsprototypen für die benötigten Funktionen,
 * zur Formulierung des SMACOF Algorithmus.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef SMACOF_HELPER_H_
#define SMACOF_HELPER_H_

#include "m.h"
#include "tm.h"
#include "utils.h"

/**
 * @brief Führt eine Guttman-Transformation auf die gewohnte Art und Weise
 * durch.
 *
 * Es wird eine Guttman-Transformation durchgeführt. Die errechnete update
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
m_t *guttmanTransformation(const tm_t *delta, const m_t *z, const tm_t *w,
                           const tm_t *pinv);

/**
 * @brief Führt eine Guttman-Transformation ohne eigenes Speichermanagement
 * durch.
 *
 * Es wird eine Guttman-Transformation durchgeführt, ohne das diese Funktion
 * ein Speichermanagement durchführt. Wird w = NULL übergeben,
 * so kann die einfachere Updateformel verwendet werden, sind Gewichte
 * angegeben, d.h. w != NULL, so wird die erweiterte Updateformel verwendet.
 *
 * @param delta Unähnlichkeitsmatrix
 * @param z aktuelle Konfiguration
 * @param x Speicherstelle an die die neue Konfiguration geschrieben wird
 * @param w Gewichtsmatrix
 * @param pinv Pseudoinverses V+
 * @param dz Distanzmatrix
 * @param bz B Matrix
 */
m_t *guttmanTransformation_nomem(const tm_t *delta, m_t *z, m_t *x,
                                 const tm_t *w, const tm_t *pinv,
                                 const tm_t *dz, tm_t *bz);

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
tm_t *computeV(const tm_t *w);

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
real_t computesigma(const m_t *x, const tm_t *delta, const tm_t *w);

/**
 * @brief Berechnet den Stresswert sigma ohne eigenes Speichermanagement.
 *
 * Es wird der Stresswert einer Konfiguration berechnet ohne das eine
 * Speichermanagement durch die Funktion durchgeführt wird.
 *
 * @param dx Distanzmatrix
 * @param delta Unähnlichkeismatrix
 * @param w Gewichtsmatrix
 * @return sigma Stress der Konfiguration
 */
real_t computesigma_nomem(const tm_t *dx, const tm_t *delta, const tm_t *w);

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
tm_t *computeBofZ(const tm_t *delta, const tm_t *dz, const tm_t *w);

/**
 * @brief Berechnet die B Matrix aus der aktuellen Distanzmatrix ohne
 * Speichermanagement.
 *
 * Es wird die B Matrix aus der aktuellen Distanzmatrix und der
 * Unähnlichkeitsmatrix berechnet. Evtl. (dev_w != NULL) erfolgt eine
 * Gewichtung. Die Funktion führt kein eigenes Speichermanagement durch.
 *
 * @param delta Unähnlichkeitsmatrix
 * @param dz Distanzmatrix
 * @param w Gewichtsmatrix
 * @param b Matrix B in die das Ergebnis geschrieben werden soll
 */
void computeBofZ_nomem(const tm_t *delta, const tm_t *dz, const tm_t *w,
                       tm_t *b);

/**
 * @brief Erzeugt eine Zufallsmatrix.
 *
 * Es wird eine Zufallsmatrix entsprechender Größe erzeugt. Bei den Einträgen
 * handelt es sich um reelle Zahlen zwischen min und max.
 *
 * @param rows Anzahl der Zeilen
 * @param cols Anzahl der Spalten
 * @param min minimaler Größe der Zahl
 * @param max maximale Größe der Zahl
 * @return Zufallsmatrix
 */
m_t *generateRandM(const unsigned int rows, const unsigned int cols,
                   const unsigned int min, const unsigned int max);

/**
 * @brief Konvertiert eine Matrix in eine Dreiecksmatrix.
 *
 * Die übergebene allgemeine Matrix wird in eine Dreiecksmatrix konvertiert.
 * Dazu wird eine Dreiecksmatrix entsprechender Größe angelegt und die
 * Einträge der unteren 'Dreieckshälfte' der übergebenen in diese kopiert.
 *
 * @param tmp_w Matrix
 * @return Dreiecksmatrix
 */
tm_t *convertMtoTM(const m_t *tmp_w);

/**
 * @brief Erzeugt eine Unähnlichkeitsmatrix mit Werten eines Beispiels.
 *
 * @return Unähnlichkeitsmatrix Beispiel
 */
tm_t *get_test_delta();

/**
 * @brief Erzeugt eine Startkonfiguration mit Werten eines Beispiels.
 *
 * @return Startkonfiguration Beispiel
 */
m_t *get_test_x();

#endif /* SMACOF_HELPER_H_ */
