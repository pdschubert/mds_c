/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file linalg.h
 * @brief Funktionsprototypen zu den mathematischen Operationen auf den
 * Matrixobjekten.
 *
 * Diese Datei enthält die Funktionsprototypen zur Realisierung der
 * mathematischen Operationen auf den allgemeinen- und Dreiecksmatrizen.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef LINALG_H_
#define LINALG_H_

#include "m.h"
#include "tm.h"
#include "utils.h"

/**
 * @brief Realisiert eine allgemeine Matrixmultiplikation.
 *
 * Die Funktion realisiert die Berechnung eines Matrixproduktes zwischen
 * zwei allgemeinen Matrizen.
 *
 * @param a allgemeine Matrix
 * @param b allgemeine Matrix
 * @return Matrixprodukt
 */
m_t *MmultM(const m_t *a, const m_t *b);

/**
 * @brief Addiert zwei allgemeine Matrizen gleicher Dimension.
 *
 * Es wird eine Addition zweier allgemeiner Matrizen gleicher Dimension
 * durchgeführt.
 *
 * @param a allgemeine Matrix
 * @param b allgemeine Matrix
 * @return Matrixsumme
 */
m_t *MaddM(const m_t *a, const m_t *b);

/**
 * @brief Berechnet ein Matrixvielfaches.
 *
 * Die Funktion berechnet das Vielfache einer allgmeinen Matrix.
 *
 * @param m allgemeine Matrix
 * @param x skalarer Wert
 * @return Matrixvielfaches
 */
m_t *Mmultscalar(const m_t *m, const real_t x);

/**
 * @brief Berechnet eine Matrixvielfaches ohne eigenes Speichermanagement.
 *
 * Die Funktion berechnet ein Matrixvielfaches. Die Multiplikation wird
 * direkt mit dem Speicher der übergebenen Matrix durchgeführt und es erfolgt
 * keine neue Speicherallokation.
 *
 * @param m allgemeine Matrix
 * @param x skalarer Wert
 */
void Mmultscalar_nomem(m_t *m, const real_t x);

/**
 * @brief Addiert einen skalaren Wert auf jeden Eintrag der Matrix.
 *
 * Die Funktion addiert den gegebenen skalaren Wert auf jeden Eintrag
 * der Matrix.
 *
 * @param m allgemeine Matrix
 * @param x skalarer Wert
 * @return Ergebnismatrix
 */
m_t *Maddscalar(const m_t *m, const real_t x);

/**
 * @brief Addiert zwei Dreiecksmatrizen gleicher Dimension.
 *
 * Es wird die Matrixsummer zweier Dreiecksmatrizen gleicher Dimension
 * gebildet.
 *
 * @param a Dreiecksmatrix
 * @param b Dreiecksmatrix
 * @return Matrixsumme
 */
tm_t *TMaddTM(const tm_t *a, const tm_t *b);

/**
 * @brief Berechnet das Vielfache einer Dreiecksmatrix.
 *
 * Es wird das Vielfache einer Dreiecksmatrix berechnet.
 *
 * @param m Dreiecksmatrix
 * @param x skalarer Wert
 */
tm_t *TMmultscalar(const tm_t *m, const real_t x);

/**
 * @brief Addiert einen skalaren Wert auf jeden Eintrag der Matrix.
 *
 * Die Funktion addiert den gegebenen skalaren Wert auf jeden Eintrag
 * der Dreiecksmatrix.
 *
 * @param m Dreiecksmatrix
 * @param x skalarer Wert
 * @return Ergebnismatrix
 */
tm_t *TMaddscalar(const tm_t *m, const real_t x);

/**
 * @brief Berechnet ein Matrixprodukt aus einer allgemeinen und einer
 * Dreiecksmatrix.
 *
 * Es wird das Matrixprodukt zwischen einer allgemeinen und einer Dreiecksmatrix
 * gebildet. Das Ergebnis wird in einer allgemeinen Matrix zurückgegeben.
 *
 * @param a allgemeine Matrix
 * @param b Dreiecksmatrix
 * @return Matrixprodukt
 */
m_t *MmultTM(const m_t *a, const tm_t *b);

/**
 * @brief Es wird das Matrixprodukt aus einer Dreiecks- und einer allgemeinen
 * Matrix berechnet.
 *
 * Es wird das Matrixprodukt as einer Dreiecksmatrix un einer allgemeinen Matrix
 * gebildet. Das Ergebnis wird in einer allgemeinen Matrix gespeichert und
 * zurückgegeben.
 *
 * @param a Dreieckmatrix
 * @param b allgemeine Matrix
 * @return Matrixprodukt
 */
m_t *TMmultM(const tm_t *a, const m_t *b);

/**
 * @brief Berechnet das Produkt aus einer Dreiecks- und einer allgemeinen Matrix
 * ohne eigenes Speichermanagement.
 *
 * Es wird das Matrixprodukt aus einer Dreiecks- und einer allgemeinen Matrix
 * berechnet. Dabei wird von dieser Funktion das Speichermanagement nicht
 * übernommmen.
 *
 * @param a Dreiecksmatrix
 * @param b allgemeine Matrix
 * @param c Adresse an die das Ergebnis geschrieben wird.
 */
void TMmultM_nomem(const tm_t *a, const m_t *b, m_t *c);

#endif /* LINALG_H_ */
