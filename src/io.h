/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file io.h
 * @brief Funktionsprototypen zur Dateiein- und ausgabe.
 *
 * Diese Datei enthält die Funktionsprototypen zur Realisierung der Dateiein-
 * und Dateiausgaben, die im Rahmen des SMACOF Algorithmus benötigt werden.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef IO_H_
#define IO_H_

#include "m.h"
#include "tm.h"

/**
 * @brief Liest eine CSV Datei ein.
 *
 * Die Funktion liest eine CSV Datei, die aus Spaltenvektoren besteht, ein und
 * interpretiert die erste Zeile der Datei als Headerzeile die in header
 * gespeichert wird. Der daraf folgende Teil wird in der Matrix data
 * gepspeichert. Der benötigte Speicher wird von der Funktion selbst ermittelt
 * und angefordert. Nach dem füllen der angeforderten Strukturen werden diese in
 * die übergebenen Adressen von data und header eingetragen. Dadurch können
 * Header- und Datenteil gleichzeitig aus der Funktion zurückgegeben werden. Die
 * Datei wird zeilenweise eingelesen.
 *
 * @warning Es wird die getline() Funktion benutzt, die den GCC Compiler
 * erfordert
 * @param path Dateipfad
 * @param data Datenmatrix
 * @param header Header
 */
void readCSVtoM(const char *path, m_t **data, char ***header);

/**
 * @brief Schreibt eine Datenmatrix mit zugehöriger Headerzeile in eine Datei.
 *
 * Die Funktion schreibt eine Datenmatrix mit zugehöriger Headerzeile in eine
 * valide CSV Datei, die die Daten mit zugehörigem Head als Spaltenvektor
 * speichert. Falls die Datei, welche durch den Dateipfad identifiziert wird,
 * bereits existiert, werden die Daten an das Dateiende angehängt. Die
 * Ergebnisse werden zeilenweise in die Datei geschrieben.
 *
 * @warning Es wird die getline() Funktion benutzt, die den GCC Compiler
 * erfordert
 * @param path Dateipfad
 * @param data Datenmatrix
 * @param header Header zur Datenmatrix
 */
void writeMtoCSV(const char *path, const m_t *data, char **header);

/**
 * @brief Liest eine CSV Datei in eine Datenmatrix.
 *
 * Liest eine CSV Datei ein, die aus einer Reihe von Spaltenvektoren besteht
 * ein. Nach dem Einlesen wird die erstellte Datenmatrix in die übergebene
 * Adresse eingetragen. Da eine Gewichtsmatrix eingelesen werden soll, muss die
 * einzulesende CSV Datei eine quadratische, symmetrische Matrix enthalten. Die
 * Ergebnisse werden zweilenweise eingelesen.
 *
 * @param path Dateipfad
 * @param data Datenmatrix
 */
void readWeightsFile(const char *path, m_t **data);

#endif /* IO_H_ */
