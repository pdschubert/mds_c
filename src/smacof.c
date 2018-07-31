/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/**
 * @file smacof.c
 * @brief Enthält die Implementation der Funktionen zur Berechnung des SMACOF
 * Algorithmus aus smacof.h.
 *
 * Die Datei enthält die Implementation des SMACOF Algorithmus in C in der
 * Standard- und in der Memory-Pool-Version.
 *
 * @author Philipp D. Schubert
 * @bug Keine Bugs bekannt.
 */

#include "smacof.h"
#include "cuda_disfuncs.cuh"
#include "cuda_linalg.cuh"
#include "cuda_smacof_helper.cuh"
#include "disfuncs.h"
#include "linalg.h"
#include "m.h"
#include "smacof_helper.h"
#include "time.h"
#include "tm.h"
#include "utils.h"
#include <math.h>
#include <string.h>

m_t *computeSMACOF(const tm_t *delta, const tm_t *w, const tm_t *pinv,
                   const real_t epsilon, const unsigned int dim,
                   const unsigned int maxiter, const unsigned int print) {
// the complete algorithm:
#if DEBUG == 1
  // make me a test delta and a test x from textbook example
  tm_t *test_delta = get_test_delta();
  m_t *test_x = get_test_x();
  printf("sigma of test-data is:%f\n", computesigma(test_x, test_delta, w));
  m_t *z = test_x;
  delta = test_delta;
#else
  // here comes the core algorithm
  m_t *z = generateRandM(delta->size, dim, 0, 10);
#endif
  m_t *x = z;
#if CUDA == 0
  real_t sigma = computesigma(x, delta, w);
#else
  real_t sigma = cuda_computesigma(x, delta, w);
#endif
  real_t sigma_prev = 0;
  real_t sigdiff;
  double iterationtime;
  unsigned int k = 0;
  while ((sigdiff = fabs(sigma_prev - sigma)) > epsilon && k < maxiter) {
    if (print)
      start_watch();
    ++k;
    sigma_prev = sigma;
    z = x;
#if CUDA == 0
    x = guttmanTransformation(delta, z, w, pinv);
#else
    x = cuda_guttmanTransformation(delta, z, w, pinv);
#endif
#if CUDA == 0
    sigma = computesigma(x, delta, w);
#else
    sigma = cuda_computesigma(x, delta, w);
#endif
#if DEBUG == 1
    if (k > 1) // for k <= 1, z does not point to dynamic memory for the test
               // case!!!
#endif
      freeM(z);
    z = NULL;
    if (print) {
      iterationtime = stop_watch();
      printf("used %f ms in iteration %u with sigma difference %f\n",
             iterationtime, k, sigdiff);
    }
  }
#if DEBUG == 1
  freeTM(test_delta);
  test_delta = NULL;
  freeM(test_x);
  test_x = NULL;
#endif
  printf("used %u iterations.\n", k);
  return x;
}

m_t *computeSMACOF_mempool(const tm_t *delta, const tm_t *w, const tm_t *pinv,
                           const real_t epsilon, const unsigned int dim,
                           const unsigned int maxiter,
                           const unsigned int print) {
  // here comes the core algorithm
  m_t *z = generateRandM(delta->size, dim, 0, 10);
  m_t *x = initM(delta->size, dim);
  memcpy(x->elems, z->elems, z->num_elems * sizeof(real_t));
  m_t *tmp = NULL;
  tm_t *d = initTM(delta->size);
  tm_t *bz = initTM(delta->size);
  calcDistanceMatrix_nomem(x, d);
  real_t sigma = computesigma_nomem(d, delta, w);
  real_t sigma_prev = 0;
  real_t sigdiff;
  double iterationtime;
  unsigned int k = 0;
  while ((sigdiff = fabs(sigma_prev - sigma)) > epsilon && k < maxiter) {
    if (print)
      start_watch();
    ++k;
    sigma_prev = sigma;
    tmp = z;
    z = x;
    x = tmp;
    x = guttmanTransformation_nomem(delta, z, x, w, pinv, d, bz);
    calcDistanceMatrix_nomem(x, d);
    sigma = computesigma_nomem(d, delta, w);
    if (print) {
      iterationtime = stop_watch();
      printf("used %f ms in iteration %u with sigma difference %f\n",
             iterationtime, k, sigdiff);
    }
  }
  freeM(z);
  freeTM(d);
  freeTM(bz);
  printf("used %u iterations.\n", k);
  return x;
}

void test_cuda_functions() {
  // m_t* x = get_test_x();
  // tm_t* delta = get_test_delta();
  // real_t sigma = computesigma(x, delta, NULL);
  // real_t cusigma = cuda_computesigma(x, delta, NULL);
  // printf("sigma:%f\n", sigma);
  // printf("cusigma:%f\n", cusigma);
  // m_t* v = initM(1000, 8);
  // for(unsigned int i = 0; i < v->num_elems; ++i) {
  //	v->elems[i] = i+1;
  //}
  // tm_t* dist = calcDeltaMatrix(v, EUCLIDEAN, 2);
  // tm_t* cudist = cuda_distanceMatrix(v);
  // printTM(dist);
  // printTM(cudist);
  // m_t* y = initM(10000000, 2);
  // for(unsigned int i = 0; i < y->num_elems; ++i)
  //	y->elems[i] = 1;
  // m_t* mult = Mmultscalar(y, 2);
  // m_t* cumult = cuda_Mmultscalar(y, 2);
  // printM(mult);
  // puts("-------");
  // printM(cumult);
  // tm_t* db = initTM(12000);
  // for(unsigned int i = 0; i < db->num_elems; ++i)
  //	db->elems[i] = i+1;
  // tm_t* b = computeBofZ(delta, dist, NULL);
  // tm_t* cub = cuda_computeBofZ(delta, dist, NULL);
  // printTM(b);
  // printTM(cub);
  // for(unsigned int i = 0; i < db->size; ++i)
  //	printf("%f\n", cub->elems[i*(i+1)/2 + i]);
  // tm_t* t = initTM(5);
  // for(unsigned int i = 0; i < t->num_elems; ++i)
  //	t->elems[i] = i + 1;
  // m_t* m = initM(5, 5);
  // for(unsigned int i = 0; i < m->num_elems; ++i)
  //	m->elems[i] = i + 1;
  // m_t* r = TMmultM(t, m);
  // m_t* cur = cuda_TMmultM(t, m);
  // printM(r);
  // puts("------------");
  // printM(cur);

  // tm_t* delta = initTM(10);
  // m_t* d = initM(1000, 2);
  // for(unsigned int i = 0; i < delta->size; ++i) {
  //	for(unsigned int j = 0; j < i; ++j)
  //	delta->elems[(i*(i+1))/2 + j] = 1;
  //}
  // for(unsigned int i = 0; i < d->num_elems; ++i) {
  //	d->elems[i] = i+1;
  //}
  // tm_t* cpuddelta = calcDeltaMatrix(d, EUCLIDEAN, 2);
  // tm_t* gpudist = cuda_distanceMatrix(d);
  // printTM(cpuddelta);
  // puts("-----------------------");
  // printTM(gpudist);
  // real_t cpusig = computesigma(d, delta, NULL);
  // real_t gpusig = cuda_computesigma(d, delta, NULL);
  // printf("cpu sigma: %f\n", cpusig);
  // printf("gpu sigma: %f\n", gpusig);
  // m_t* cpumult = TMmultM(delta, d);
  // m_t* gpumult = cuda_TMmultM(delta, d);
  // printM(cpumult);
  // printM(gpumult);
}
