/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/** @file cuda_smacof.cu
 *  @brief Implementation der Prototypen aus cuda_smacof.cuh.
 *
 *  Hier ist die Wrapperfunktion zur Berechnung des SMACOF Algorithmus
 *  in der Memory Pool Version implementiert.
 *
 *  @author Philipp D. Schubert
 *  @bug Keine Bugs bekannt.
 */

#include <math.h>
#include <stdlib.h>

extern "C" {
#include "cuda_disfuncs.cuh"
#include "cuda_linalg.cuh"
#include "cuda_smacof_helper.cuh"
#include "m.h"
#include "smacof_helper.h"
#include "time.h"
#include "tm.h"
#include "utils.cuh"
#include "utils.h"
}

extern "C" m_t *cuda_computeSMACOF_mempool(const tm_t *delta, const tm_t *w,
                                           const tm_t *pinv,
                                           const real_t epsilon,
                                           const unsigned int dim,
                                           const unsigned int maxiter,
                                           const unsigned int print) {
  // here comes the core algorithm
  m_t *z = generateRandM(delta->size, dim, 0, 10);
  m_t *x = initM(delta->size, dim);
  memcpy(x->elems, z->elems, z->num_elems * sizeof(real_t));
  real_t *dev_z, *dev_x, *dev_tmp = NULL, *dev_delta, *dev_pinv, *dev_w = NULL,
                         *dev_d, *dev_bz;
  CUERROR(cudaMalloc((void **)&dev_z, z->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_z, z->elems, z->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  CUERROR(cudaMalloc((void **)&dev_x, x->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_x, x->elems, x->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  CUERROR(cudaMalloc((void **)&dev_delta, delta->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_delta, delta->elems, delta->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  if (w != NULL) {
    CUERROR(cudaMalloc((void **)&dev_w, w->num_elems * sizeof(real_t)));
    CUERROR(cudaMemcpy(dev_w, w->elems, w->num_elems * sizeof(real_t),
                       cudaMemcpyHostToDevice));
    CUERROR(cudaMalloc((void **)&dev_pinv, pinv->num_elems * sizeof(real_t)));
    CUERROR(cudaMemcpy(dev_pinv, pinv->elems, pinv->num_elems * sizeof(real_t),
                       cudaMemcpyHostToDevice));
  }
  CUERROR(cudaMalloc((void **)&dev_d, delta->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_bz, delta->num_elems * sizeof(real_t)));

  cuda_distanceMatrix_nomem(dev_x, x->rows, x->cols, dev_d);
  real_t sigma =
      cuda_computesigma_nomem(dev_d, delta->num_elems, dev_delta, dev_w);
  real_t sigma_prev = 0;
  real_t sigdiff;
  double iterationtime;
  unsigned int k = 0;
  while ((sigdiff = fabs(sigma_prev - sigma)) > epsilon && k < maxiter) {
    if (print)
      start_watch();
    ++k;
    sigma_prev = sigma;
    dev_tmp = dev_z;
    dev_z = dev_x;
    dev_x = dev_tmp;
    dev_x = cuda_guttmanTransformation_nomem(
        dev_delta, delta->size, delta->num_elems, dev_z, z->rows, z->cols,
        dev_x, dev_w, dev_pinv, dev_d, dev_bz);
    cuda_distanceMatrix_nomem(dev_x, x->rows, x->cols, dev_d);
    sigma = cuda_computesigma_nomem(dev_d, delta->num_elems, dev_delta, dev_w);
    if (print) {
      iterationtime = stop_watch();
      printf("used %f ms in iteration %u with sigma difference %f\n",
             iterationtime, k, sigdiff);
    }
  }
  freeM(z);
  cudaFree(dev_d);
  cudaFree(dev_bz);
  cudaFree(dev_delta);
  if (w != NULL) {
    cudaFree(dev_w);
    cudaFree(dev_pinv);
  }
  CUERROR(cudaMemcpy(x->elems, dev_x, x->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_x);
  cudaFree(dev_z);
  printf("used %u iterations.\n", k);
  return x;
}
