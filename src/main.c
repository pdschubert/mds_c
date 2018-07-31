/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/
// ============================================================================
//  Description : Compute the SMACOF algorithm in OpenMP and CUDA C
//  ===========================================================================

#include "cuda_disfuncs.cuh"
#include "cuda_linalg.cuh"
#include "cuda_smacof.cuh"
#include "cuda_smacof_helper.cuh"
#include "disfuncs.h"
#include "io.h"
#include "linalg.h"
#include "m.h"
#include "pseudoinverse.hh"
#include "smacof.h"
#include "smacof_helper.h"
#include "tm.h"
#include "utils.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  // parameters for this programm
  char *inputfilepath;
  char *outputfilepath;
  unsigned int disfunction;
  unsigned int maxiteration;
  real_t epsilon;
  unsigned int resultdim;
  unsigned int minkowski_p;
  unsigned int print;
  char *weightsfilepath = NULL;
  // fill parameters with test values
#if DEBUG == 1
  inputfilepath = "/home/philipp/GIT-Repos/MultidimensionalScaling/mds_c/"
                  "example-data/example_data.csv";
  outputfilepath = "/home/philipp/GIT-Repos/MultidimensionalScaling/mds_c/"
                   "build/results-example_data.csv";
  disfunction = EUCLIDEAN;
  maxiteration = 1000;
  epsilon = 0.0;
  resultdim = 2;
  minkowski_p = 2;
  print = 1;
  // weightsfilepath = "/home/philipp/Schreibtisch/myweights.txt";
  // fill parameters with user values
#elif DEBUG == 0
  const char *usage =
      "####################################################################\n"
      "# SMACOF - Philipp D. Schubert - philipp@it-schubert.com           #\n"
      "####################################################################\n"
      "usage: <prog> <input> <output> <disfunc> <maxiter> <epsilon> "
      "<resultdim> <metric_p> <print>\n"
      "parameter explanation:\n"
      "\tinput - a valid path to csv file to analyse: path\n"
      "\touput - a filename for the outputfile: path\n"
      "\tdisfunc - dissimilarity measure: ...\n"
      "\t\t0 - euclidean\n"
      "\t\t1 - cityblock\n"
      "\t\t2 - minkowski\n"
      "\t\t3 - correlation\n"
      "\t\t4 - angularseperatin\n"
      "\t\t5 - wavehedges\n"
      "\t\t6 - bahattacharyya\n"
      "\t\t7 - soergel\n"
      "\t\t8 - braycurtis\n"
      "\t\t9 - divergence\n"
      "\t\t10 - canberra\n"
      "\t\t11 - none\n"
      "\tmaxiter - maximum number of iterations to use: a positive integer\n"
      "\tepsilon - the maximum error value which is allowed: a positive real_t\n"
      "\tresultdim - number of dimensions for result vectors: a positive "
      "integer\n"
      "\tminkowski_p - a value to adjust the minkowski distance: a positive "
      "integer\n"
      "\tprint - be verbose: 0 or 1\n"
      "\tweights - (optinal) a valid path to csv file containing the weights\n";
  if (argc < 9) {
    printf("ERROR: wrong number of arguments\n");
    printf("%s\n", usage);
    exit(-1);
  }
  puts("---------- program arguments ----------");
  inputfilepath = argv[1];
  printf("inputfilepath: %s\n", inputfilepath);
  outputfilepath = argv[2];
  printf("outputfilepath: %s\n", outputfilepath);
  disfunction = atoi(argv[3]);
  printf("dissimilarity function: %u\n", disfunction);
  maxiteration = atoi(argv[4]);
  printf("max iteration: %u\n", maxiteration);
  epsilon = atof(argv[5]);
  printf("epsilon: %f\n", epsilon);
  resultdim = atoi(argv[6]);
  printf("resultdim: %u\n", resultdim);
  minkowski_p = atoi(argv[7]);
  printf("metric p: %u\n", minkowski_p);
  print = atoi(argv[8]);
  printf("print: %u\n", print);
  if (argc == 10) {
    weightsfilepath = argv[9];
    printf("weightsfilepath: %s\n", argv[9]);
  }
  if (epsilon < 0.0) {
    printf("epsilon has to be non-negativ\n");
    printf("%s\n", usage);
    exit(-1);
  }
#endif

  m_t *data;
  tm_t *weights = NULL;
  tm_t *pinv = NULL;
  char **header;
  puts("read input file ...");
  readCSVtoM(inputfilepath, &data, &header);

  puts("compute delta matrix ...");
  tm_t *delta = calcDeltaMatrix(data, disfunction, minkowski_p);
  freeM(data);
  data = NULL;

  if (weightsfilepath != NULL) {
    // build up weight matrix
    puts("read weights file ...");
    m_t *tmp;
    readWeightsFile(weightsfilepath, &tmp);
    weights = convertMtoTM(tmp);
    freeM(tmp);
    tmp = NULL;
    CERROR(weights->size != delta->size,
           "dimensions of weight matrix and delta matrix does not match");
    // build up V matrix and pseudoinverse of V
    puts("compute v matrix ...");
#if CUDA == 0
    tm_t *v = computeV(weights);
#else
    tm_t *v = cuda_computeV(weights);
#endif
    puts("compute pseudoinverse of v ...");
    pinv = pseudoinverse(v);
    freeTM(v);
    v = NULL;
  }
  puts("compute the smacof solution - start iterating ...");
#if MPOOL == 1
#if CUDA == 0
  m_t *result = computeSMACOF_mempool(delta, weights, pinv, epsilon, resultdim,
                                      maxiteration, print);
#else
  m_t *result = cuda_computeSMACOF_mempool(delta, weights, pinv, epsilon,
                                           resultdim, maxiteration, print);
#endif
#else
  m_t *result = computeSMACOF(delta, weights, pinv, epsilon, resultdim,
                              maxiteration, print);
#endif
  puts("found result configuration.");
  if (print) {
    puts("result:");
    printM(result);
  }
  puts("write output file ...");
  writeMtoCSV(outputfilepath, result, header);

  // free all the stuff we have used
  puts("free used memory ...");
  for (unsigned int i = 0; i < delta->size; ++i)
    free(header[i]);
  free(header);
  header = NULL;
  freeTM(delta);
  delta = NULL;
  if (weights != NULL) {
    freeTM(weights);
    weights = NULL;
    freeTM(pinv);
    pinv = NULL;
  }
  freeM(result);
  result = NULL;

  // testing cuda functions:
  // puts("--- CUDA TESTS -----------------");
  // test_cuda_functions();

  puts("program finished.");
  printf("check outputfile: %s for results\n", outputfilepath);
  exit(0);
}
