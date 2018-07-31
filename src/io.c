/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/** @file io.c
 *  @brief Implementation der Prototypen aus io.h.
 *
 *  In dieser Datei sind die Funktionen zur Dateiein- und ausgabe, die zur
 *  Durchführung des SMACOF Algorithmus benötigt werden.
 *
 *  @author Philipp D. Schubert
 *  @bug Keine Bugs bekannt.
 */

#define _GNU_SOURCE
#include "io.h"
#include "m.h"
#include "tm.h"
#include "utils.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void readCSVtoM(const char *path, m_t **data, char ***header) {
  /* caution: getline() function only available for GNU gcc
   * getline() does the memory allocation for our buffers automatically,
   * if all buffers given to the function are initialized with ZERO!
   * see getline() documentation for details, it is one of the
   * most secure file-reading functions and should be used.
   */
  FILE *file;
  if ((file = fopen(path, "r")) == NULL) {
    perror("could not open file for reading data");
    HEREANDNOW;
    exit(-1);
  }

  // count the number of lines and number of attributes in file
  char buf;
  int number_lines = 1;
  int number_attributes = 0;
  for (;;) {
    size_t n = fread(&buf, 1, 1, file);
    if (n < 1)
      break;
    if (buf == '\n')
      ++number_lines;
    if (buf == ',')
      ++number_attributes;
  }
  number_attributes /= number_lines;
  number_attributes++;
  printf("number of attributes in csv: %d\n", number_attributes);
  printf("number of lines in csv: %d\n", number_lines);
  // set file descriptor to beginning of file
  rewind(file);
  // read the first line and interpret it as the csv header
  *header = (char **)malloc(number_attributes * sizeof(char *));
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  const char delimiter[] = ",;";
  char *tok;
  int k = 0;
  if ((read = getline(&line, &len, file)) != -1) {
    tok = strtok(line, delimiter);
    while (tok != NULL) {
      (*header)[k] = strdup(tok);
      tok = strtok(NULL, delimiter);
      ++k;
    }
  }
  // setup our matrix where to read the data in
  // 'number_lines - 1' because we previously read the header line
  *data = initM(number_attributes, number_lines - 1);
  int i = 0;
  int j = 0;
  while ((read = getline(&line, &len, file)) != -1) {
    tok = strtok(line, delimiter);
    while (tok != NULL) {
      (*data)->elems[j * (*data)->cols + i] = strtof(tok, NULL);
      tok = strtok(NULL, delimiter);
      ++j;
    }
    j = 0;
    ++i;
  }
  CERROR(0 != fclose(file), "could not close file properly");
  free(line);
}

void writeMtoCSV(const char *path, const m_t *data, char **header) {
  FILE *file = fopen(path, "a+");
  CERROR(file == NULL, "could not open file for writing results");

  // write header first
  for (unsigned int i = 0; i < data->rows; ++i) {
    fputs(header[i], file);
    if (i < data->rows - 1)
      fputs(",", file);
  }
  // here comes the data
  char buffer[16];
  for (unsigned int i = 0; i < data->cols; ++i) {
    for (unsigned int j = 0; j < data->rows; ++j) {
      // printf("indices - i:%d, j:%d value:%f\n", j, i, data->elems[j *
      // data->cols + i]);
      snprintf(buffer, 16, "%f", data->elems[j * data->cols + i]);
      fputs(buffer, file);
      if (j < data->rows - 1)
        fputs(",", file);
      else
        fputs("\n", file);
    }
  }
  CERROR(0 != fclose(file), "could not close file properly");
}

void readWeightsFile(const char *path, m_t **weights) {
  /* caution: getline() function only available for GNU gcc
   * getline() does the memory allocation for our buffers automatically,
   * if all buffers given to the function are initialized with ZERO!
   * see getline() documentation for details, it is one of the
   * most secure file-reading functions and should be used.
   */
  FILE *file;
  int nret;
  char **line = malloc(sizeof(char *));
  *line = NULL;
  size_t *t = malloc(0);

  if ((file = fopen(path, "r")) == NULL) {
    perror("could not open file for reading weights");
    HEREANDNOW;
    exit(-1);
  }

  // count the number of lines and number of attributes in file
  char ch;
  int number_lines = 0;
  int number_attributes = 0;
  while (!feof(file)) {
    ch = fgetc(file);
    if (ch == '\n') {
      number_lines++;
    }
    if (ch == ',') {
      number_attributes++;
    }
  }
  number_attributes /= number_lines;
  number_attributes++;

  // set file descriptor to beginning of file
  rewind(file);

  char delimiter[] = ",;";
  // setup our matrix where to read the data in
  *weights = initM(number_attributes, number_lines);
  int i = 0;
  int j = 0;
  char *tok;
  while ((nret = getline(line, t, file)) > 0) {
    tok = strtok(*line, delimiter);
    while (tok != NULL) {
      (*weights)->elems[j * (*weights)->cols + i] = strtof(tok, NULL);
      tok = strtok(NULL, delimiter);
      ++j;
    }
    j = 0;
    ++i;
  }

  CERROR(0 != fclose(file), "could not close file properly");

  free(*line);
  free(line);
  free(t);
}
