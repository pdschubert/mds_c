/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file disfuncs.h
 * @brief Funktionsprototypen zur Berechnung der Unähnlichkeitsmatrix, sowie der
 * Distanzmatrix.
 *
 * Diese Datei enthält die Funktionsprototypen zur Berechnun der
 * Unähnlichkeitsmatrix bzw. Distanzmatrix in der 'nomem' Version. Desweiteren
 * sind die einzelnen Unähnlichkeitsfunktionen hier als inline-Funktionen
 * definiert, da diese in der Regel sehr klein sind und sehr häufig benötigt
 * werden.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef DISFUNCS_H_
#define DISFUNCS_H_

#include "m.h"
#include "tm.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Berechnet aus einer Datenmatrix die Unähnlichkeitsmatrix.
 *
 * Diese Funktion berechnet aus einer gegebenen Datenmatrix die
 * quadratische Unähnlichkeitsmatrix Delta. Dabei wird die gewählte
 * Unähnlichkeitsfunktion, sowie der Parameter p verwendet. P wird für
 * einige Unähnlichkeitsfunktionen benötigt.
 *
 * @param data Datenmatrix
 * @param d Unähnlichkeitsfunktion
 * @param p Parameter für Unähnlichkeitsfunktionen.
 */
tm_t *calcDeltaMatrix(const m_t *data, const enum delta d, const int p);

/**
 * @brief Berechnet die (euklidische) Distanzmatrix aus gegebener Datenmatrix.
 *
 * Diese Funktion berechnet aus der Datenmatrix die euklidische Distanzmatrix.
 * Das Speichermanagement wird von dieser Funktion nicht selbstständig
 * übernommen.
 *
 * @param x Datenmatrix
 * @param dx Matrix, in die die Distanzmatrix geschrieben werden soll.
 */
void calcDistanceMatrix_nomem(const m_t *x, tm_t *dx);

/**
 * @brief Berechnet das Maximum zweier Gleitkommazahlen.
 *
 * @param a Gleitkommazahl
 * @param b Gleitkommazahl
 * @return Maximum
 */
static inline real_t mymax(const real_t a, const real_t b) {
  return (a > b) ? a : b;
}

/**
 * @brief Berechnet das Minimum zweier Gleitkommazahlen.
 *
 * @param a Gleitkommazahl
 * @param b Gleitkommazahl
 * @return Minimum
 */
static inline real_t mymin(const real_t a, const real_t b) {
  return (a < b) ? a : b;
}

/**
 * @brief Berechnet den Mittelwert eines Arrays.
 *
 * @param a Array
 * @param len Länge des Arrays
 * @return Mittelwert
 */
static inline real_t mean(const real_t *a, const unsigned int len) {
  real_t mean = 0;
  for (unsigned int i = 0; i < len; i++) {
    mean += a[i];
  }
  return mean / (real_t)len;
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der CANBERRA Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t canberra(const real_t *__restrict a,
                              const real_t *__restrict b,
                              const unsigned int len) {
  real_t sum = 0;
  for (unsigned int i = 0; i < len; ++i) {
    if (a[i] != 0 && b[i] != 0)
      sum += fabsf(a[i] - b[i]) / (fabsf(a[i]) + fabsf(b[i]));
  }
  return sum;
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der DIVERGENCE Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t divergence(const real_t *__restrict a,
                                const real_t *__restrict b,
                                const unsigned int len) {
  real_t sum = 0;
  real_t den = 0;
  for (unsigned int i = 0; i < len; ++i) {
    den = (a[i] + b[i]) * (a[i] + b[i]);
    sum += (den != 0.0) ? ((a[i] - b[i]) * (a[i] - b[i])) / den : 0.0;
  }
  return sum;
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der BRAYCURTIS Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t braycurtis(const real_t *__restrict a,
                                const real_t *__restrict b,
                                const unsigned int len) {
  real_t num = 0;
  real_t den = 2;
  for (unsigned int i = 0; i < len; ++i) {
    num += fabsf(a[i] - b[i]);
    den += a[i] + b[i];
  }
  return num / den;
}

/**
 * @brief Berechent einen Unähnlichkeitswert entsprechend der SOERGEL Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t soergel(const real_t *__restrict a,
                             const real_t *__restrict b,
                             const unsigned int len) {
  real_t num = 0;
  real_t den = 0;
  for (unsigned int i = 0; i < len; ++i) {
    num += fabsf(a[i] - b[i]);
    den += mymax(a[i], b[i]);
  }
  return (den != 0.0) ? num / den : 0.0;
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der BAHATTACHARYYA
 * Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vekotrs
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t bahattacharyya(const real_t *__restrict a,
                                    const real_t *__restrict b,
                                    const unsigned int len) {
  real_t sum = 0;
  for (unsigned int i = 0; i < len; ++i) {
    CERROR(a[i] < 0.0 || b[i] < 0.0,
           "Negative numbers for bahattacharyya distance not allowed!");
    sum += (sqrtf(a[i]) - sqrtf(b[i])) * (sqrtf(a[i]) - sqrtf(b[i]));
  }
  return sqrtf(sum);
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der WAVEHEDGES Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t wavehedges(const real_t *__restrict a,
                                const real_t *__restrict b,
                                const unsigned int len) {
  real_t sum = 0;
  for (unsigned int i = 0; i < len; ++i) {
    sum += (mymax(a[i], b[i]) != 0.0)
               ? 1 - (mymin(a[i], b[i]) / mymax(a[i], b[i]))
               : 1.0;
  }
  return sum;
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der ANGULARSPERATION
 * Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vekotrs
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t angularseperation(const real_t *__restrict a,
                                       const real_t *__restrict b,
                                       const unsigned int len) {
  real_t num = 0;
  real_t den_a = 0;
  real_t den_b = 0;
  real_t den = 0;
  for (unsigned int i = 0; i < len; ++i) {
    num += a[i] * b[i];
    den_a += a[i] * a[i];
    den_b += b[i] * b[i];
  }
  den = sqrtf(den_a * den_b);
  return (den != 0.0) ? 1.0 - (num / den) : 1.0;
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der CORRELATION
 * Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t correlation(const real_t *__restrict a,
                                 const real_t *__restrict b,
                                 const unsigned int len) {
  real_t num = 0;
  real_t den_a = 0;
  real_t den_b = 0;
  real_t a_mean = mean(a, len);
  real_t b_mean = mean(b, len);
  real_t den = 0;
  for (unsigned int i = 0; i < len; ++i) {
    num += (a[i] - a_mean) * (b[i] - b_mean);
    den_a += (a[i] - a_mean) * (a[i] - a_mean);
    den_b += (b[i] - b_mean) * (b[i] - b_mean);
  }
  den = sqrtf(den_a * den_b);
  return (den != 0.0) ? 1.0 - (num / den) : 1.0;
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der CITYBLOCK Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t cityblock(const real_t *__restrict a,
                               const real_t *__restrict b,
                               const unsigned int len) {
  real_t sum = 0;
  for (unsigned int i = 0; i < len; ++i) {
    sum += fabs(a[i] - b[i]);
  }
  return sum;
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der EUCLIDEAN Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t euclidean(const real_t *__restrict a,
                               const real_t *__restrict b,
                               const unsigned int len) {
  real_t sum = 0;
  for (unsigned int i = 0; i < len; ++i) {
    sum += (a[i] - b[i]) * (a[i] - b[i]);
  }
  return sqrtf(sum);
}

/**
 * @brief Berechnet einen Unähnlichkeitswert entsprechend der MINKOWSKI Formel.
 *
 * @warning Vektoren dürfen nicht überlappen
 * @param a Anfang eines Vektors
 * @param b Anfang eines Vektors
 * @param p Parameter der Unähnlichkeitsformel
 * @param len berücksichtigte Länge
 * @return Unähnlichkeitswert
 */
static inline real_t minkowski(const real_t *__restrict a,
                               const real_t *__restrict b, const int p,
                               const unsigned int len) {
  real_t sum = 0;
  for (unsigned int i = 0; i < len; ++i) {
    sum += powf(fabsf(a[i] - b[i]), (real_t)p);
  }
  return powf(sum, (real_t)1.0 / p);
}

#endif /* DISFUNCS_H_ */
