/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file time.h
 * @brief Enthält die Datenstrukturen und Funktionen, um eine Stoppuhr zu
 * realisieren.
 *
 * Diese Datei enthält die Datenstrukturen, sowie die Funktionen, die nötig
 * sind, um eine einfach Stopuhr zu realisieren. Um den Mehraufwand der
 * Funktionsaufrufe zu umgehen sind diese Funktionen inline Implementiert.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#ifndef TIME_H_
#define TIME_H_

#include <sys/time.h>

/**
 * Variablen, die die Start- und Endzeiten einer Zeitmessung speichern.
 */
static struct timeval starttime, endtime;

/**
 * @brief Ermittelt die aktuelle Systemzeit.
 *
 * @return aktuelle Systemzeit
 */
static inline struct timeval current_time() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv;
}

/**
 * @brief Berechnet mittels Start- und Endzeit die verstrichene Zeit.
 *
 * @see start_watch()
 * @see stop_watch()
 * @return Zeit die zwischen Start und Ende der Messung verstrichen ist
 */
static inline double elapsed_ms() {
  return ((endtime.tv_sec - starttime.tv_sec) * 1000.0 +
          (endtime.tv_usec - starttime.tv_usec) * 0.001);
}

/**
 * @brief Ermittelt die aktuelle Systemzeit und speichert diese in starttime.
 */
static inline void start_watch() { starttime = current_time(); }

/**
 * @brief Ermittelt die aktuelle Systemzeit und speichert diese in endtime,
 * anschließend wird die verstrichene Zeit berechnet.
 *
 * @return zwischen start und stop verstrichene Zeit
 */
static inline double stop_watch() {
  endtime = current_time();
  return elapsed_ms();
}

#endif /* TIME_H_ */
