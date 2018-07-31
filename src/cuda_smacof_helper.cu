/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/** @file cuda_smacof_helper.cu
 *  @brief Implementation der Prototypen aus cuda_smacof_helper.cuh.
 *
 *  Hier sind die Wrapperfunktionen der Prototypen aus cuda_smacof_helper.cuh
 *  implementiert, sowie entsprechende CUDA Kernel, die die Berechnungen auf der
 *  CUDA Karte durchführt.
 *
 *  @author Philipp D. Schubert
 *  @bug Keine Bugs bekannt.
 */

#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>

// my c includes
extern "C" {
#include "cuda_disfuncs.cuh"
#include "cuda_linalg.cuh"
#include "cuda_smacof_helper.cuh"
#include "m.h"
#include "tm.h"
#include "utils.cuh"
#include "utils.h"
}

static const int THREADS_PER_BLOCK = 512;
/**
 * @brief Zu benutzende Größe der TILES von den Kernel Funktionen in diesem
 * Modul.
 */
#define TILE_WIDTH 32

/**
 * @brief Test Kernel zur Überprüfung des Moduls cuda_smacof_helper.cu.
 */
__global__ void kernel_smacof_helper_test() {
  printf("cuda_smacof_helper: working fine!\n");
}

extern "C" void cuda_smacof_helper_test() {
  kernel_smacof_helper_test<<<1, 1>>>();
  cudaDeviceSynchronize();
}

extern "C" m_t *cuda_guttmanTransformation(const tm_t *delta, const m_t *z,
                                           const tm_t *w, const tm_t *pinv) {
  m_t *update = initM(z->rows, z->cols);
  const real_t norm = 1.0 / delta->size;
  // CUDA version of guttman-transformation
  // if all weights are equal to 1, use simple guttman formula
  // formula in 'linear' math: X_update = n^-1 * B(Z) * Z
  tm_t *dz = cuda_distanceMatrix(z);
  tm_t *bz = cuda_computeBofZ(delta, dz, w);
  m_t *bztimesz = cuda_TMmultM(bz, z);
  if (w == NULL) {
    update = cuda_Mmultscalar(bztimesz, norm);
  } else {
    update = cuda_TMmultM(pinv, bztimesz);
  }
  // cleanup temporary stuff
  freeTM(dz);
  dz = NULL;
  freeTM(bz);
  bz = NULL;
  freeM(bztimesz);
  bztimesz = NULL;
  return update;
}

extern "C" real_t *cuda_guttmanTransformation_nomem(
    const real_t *dev_delta, unsigned int delta_size,
    unsigned int delta_num_elems, real_t *dev_z, unsigned int zrows,
    unsigned int zcols, real_t *dev_x, const real_t *dev_w,
    const real_t *dev_pinv, const real_t *dev_d, real_t *dev_bz) {
  const real_t norm = 1.0 / delta_size;
  // CUDA version of guttman-transformation
  // if all weights are equal to 1, use simple guttman formula
  // formula in 'linear' math: X_update = n^-1 * B(Z) * Z
  cuda_computeBofZ_nomem(dev_delta, delta_size, delta_num_elems, dev_d, dev_w,
                         dev_bz);
  cuda_TMmultM_nomem(dev_bz, delta_size, dev_z, zrows, zcols, dev_x);
  if (dev_w == NULL) {
    cuda_Mmultscalar_nomem(dev_x, norm, zrows * zcols);
    return dev_x;
  } else {
    // formula in 'linear' math: X_update = V+ * B(Z) * Z
    cuda_TMmultM_nomem(dev_pinv, delta_size, dev_x, zrows, zcols, dev_z);
    return dev_z;
  }
}

/**
 * @brief Berechnet die Spalten- bzw. Zeilensummen einer Dreiecksmatrix.
 *
 * Es werden die Spalten- bzw. Zeilensummen einer Dreiecksmatrix berechnet,
 * welche anschließend negiert und in der Hauptdiagonalen der Eingabematrix
 * eingetragen wird.
 *
 * @warning GPU Kernel Funktion
 * @param v Eingabematrix
 * @param len Anzahl der Spalten von v
 */
__global__ void kernel_computeV_sum(real_t *v, const unsigned int len) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len) {
    real_t sum = 0.0;
    for (unsigned int i = 0; i < len; ++i) {
      sum +=
          (i < tid) ? v[(tid * (tid + 1)) / 2 + i] : v[(i * (i + 1)) / 2 + tid];
    }
    v[(tid * (tid + 1)) / 2 + tid] = -sum;
  }
}

/**
 * @brief Kehrt das Vorzeichen jedes Eintrags einer Matrix um.
 *
 * @warning GPU Kernel Funktion
 * @param v Eingabematrix
 * @param len Anzahl der Matrixeinträge
 */
__global__ void kernel_computeV_flip_sign(real_t *v, unsigned int len) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len)
    v[tid] = -v[tid];
}

extern "C" tm_t *cuda_computeV(const tm_t *w) {
  tm_t *v = initTM(w->size);
  real_t *dev_v;
  CUERROR(cudaMalloc((void **)&dev_v, v->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_v, w->elems, w->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  unsigned int gridsize = (v->num_elems % THREADS_PER_BLOCK == 0)
                              ? v->num_elems / THREADS_PER_BLOCK
                              : v->num_elems / THREADS_PER_BLOCK + 1;
  kernel_computeV_flip_sign<<<gridsize, THREADS_PER_BLOCK>>>(dev_v,
                                                             v->num_elems);
  kernel_computeV_sum<<<v->size, THREADS_PER_BLOCK>>>(dev_v, v->size);
  CUERROR(cudaMemcpy(v->elems, dev_v, v->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_v);
  return v;
}

/**
 * @brief Berechnet den Stresswert sigma zwischen Distanzmatrix der aktuellen
 * Konfiguration und Unähnlichkeitsmatrix.
 *
 * Diese Funktion berechnet den (gewichteten) Stresswert sigma zwischen der
 * Distanzmatrix einer Konfiguration und der Unähnlichkeitsmatrix. Dazu werden
 * eintragsweise (gewichtete) quadratische Distanzen berechnet, die anschließend
 * aufsummiert werden müssen. Die Summation der Teilwerte wir hier logarithmisch
 * (log_2) durchgeführt. Am Ende der Funktion verbleibt jeweil ein Wert im
 * geteilten Speicher cache, falls also mehrere Blöcke gestartet werden müssen,
 * so müssen die einzelnen verbleibenden Werte in cache auf der CPU (in der
 * zugehörigen Wrapperfunktion) aufaddiert werden.
 *
 * @warning GPU Kernel Funktion
 * @param dx Distanzfunktion einer Konfiguration
 * @param delta Unähnlichkeitsfunktion
 * @param w Gewichtsmatrix
 * @param partial_sigma Teilsigma, welches von einem Block berechnet werden soll
 * @param len Länge des Teilsigmas
 */
__global__ void kernel_computesigma(const real_t *dx, const real_t *delta,
                                    const real_t *w, real_t *partial_sigma,
                                    const int len) {
  __shared__ real_t cache[THREADS_PER_BLOCK];
  unsigned int bdx = blockDim.x;
  unsigned int tid = bdx * blockIdx.x + threadIdx.x;
  unsigned int cacheid = threadIdx.x;
  real_t tmp = 0;
  while (tid < len) {
    tmp += (w == NULL)
               ? (dx[tid] - delta[tid]) * (dx[tid] - delta[tid])
               : w[tid] * ((dx[tid] - delta[tid]) * (dx[tid] - delta[tid]));
    tid += gridDim.x * bdx;
  }
  cache[cacheid] = tmp;
  __syncthreads();
  unsigned int i = bdx / 2;
  while (i != 0) {
    if (cacheid < i)
      cache[cacheid] += cache[cacheid + i];
    __syncthreads();
    i /= 2;
  }
  if (cacheid == 0) {
    partial_sigma[blockIdx.x] = cache[0];
  }
}

extern "C" real_t cuda_computesigma(const m_t *x, const tm_t *delta,
                                    const tm_t *w) {
  tm_t *dx = cuda_distanceMatrix(x);
  const unsigned int gridsize =
      min(32, (dx->num_elems + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  real_t *partial_sigma = (real_t *)malloc(gridsize * sizeof(real_t));
  real_t *dev_dx, *dev_delta, *dev_w = NULL, *dev_partial_sigma;
  CUERROR(cudaMalloc((void **)&dev_dx, dx->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_delta, delta->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_partial_sigma, gridsize * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_dx, dx->elems, dx->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  CUERROR(cudaMemcpy(dev_delta, delta->elems, delta->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  if (w != NULL) {
    CUERROR(cudaMalloc((void **)&dev_w, w->num_elems * sizeof(real_t)));
    CUERROR(cudaMemcpy(dev_w, w->elems, w->num_elems * sizeof(real_t),
                       cudaMemcpyHostToDevice));
  }
  kernel_computesigma<<<gridsize, THREADS_PER_BLOCK>>>(
      dev_dx, dev_delta, dev_w, dev_partial_sigma, dx->num_elems);
  CUERROR(cudaMemcpy(partial_sigma, dev_partial_sigma,
                     gridsize * sizeof(real_t), cudaMemcpyDeviceToHost));
  real_t sigma = 0;
  for (unsigned int i = 0; i < gridsize; ++i) {
    sigma += partial_sigma[i];
  }
  freeTM(dx);
  cudaFree(dev_dx);
  cudaFree(dev_delta);
  cudaFree(dev_partial_sigma);
  free(partial_sigma);
  if (w != NULL) {
    cudaFree(dev_w);
  }
  return sigma;
}

extern "C" real_t cuda_computesigma_nomem(const real_t *dev_d,
                                          const unsigned int dnum_elems,
                                          const real_t *dev_delta,
                                          const real_t *dev_w) {
  const unsigned int gridsize =
      min(32, (dnum_elems + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  real_t *partial_sigma = (real_t *)malloc(gridsize * sizeof(real_t));
  real_t *dev_partial_sigma;
  CUERROR(cudaMalloc((void **)&dev_partial_sigma, gridsize * sizeof(real_t)));
  kernel_computesigma<<<gridsize, THREADS_PER_BLOCK>>>(
      dev_d, dev_delta, dev_w, dev_partial_sigma, dnum_elems);
  CUERROR(cudaMemcpy(partial_sigma, dev_partial_sigma,
                     gridsize * sizeof(real_t), cudaMemcpyDeviceToHost));
  real_t sigma = 0;
  for (unsigned int i = 0; i < gridsize; ++i) {
    sigma += partial_sigma[i];
  }
  cudaFree(dev_partial_sigma);
  free(partial_sigma);
  return sigma;
}

/**
 * @brief Berechnet die Quotienten der Matrix B.
 *
 * Die Funktion berechnet die (gewichteten) Quotienten aus der
 * Unähnlichkeitsmatrix und der Distanzmatrix einer Konfiguration. Fall sich der
 * Nenner zu 0.0 ergibt, so wird der gesamte Quotient als 0.0 gewertet.
 *
 * @warning GPU Kernel Funktion
 * @param delta Unähnlichkeitsmatrix
 * @param dz Distanzmatrix einer Konfiguration
 * @param w Gewichtsmatrix
 * @param b Ergebnismatrix
 * @param len Anzahl der Matrixeinträge von delta, dz und b
 * @param bsize Anzahl der Zeilen bzw. Spalten von delta, dz und b
 */
__global__ void kernel_computeB_frac(const real_t *delta, const real_t *dz,
                                     const real_t *w, real_t *b,
                                     const unsigned int len,
                                     const unsigned int bsize) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len) {
    if (w == NULL) {
      b[tid] = (dz[tid] != 0.0) ? -(delta[tid] / dz[tid]) : 0.0;
    } else {
      b[tid] = (dz[tid] != 0.0) ? -((w[tid] * delta[tid]) / dz[tid]) : 0.0;
    }
  }
}

/**
 * @brief Berechnet die Spalten bzw. Zeilensumme einer Dreiecksmatrix.
 *
 * Es wird die Spalten bzw. Zeilensumme einer Dreiecksmatrix berechnet, die
 * anschließend negiert in die Hauptdiagonal der Matrix eingetragen wird.
 *
 * @warning GPU Kernel Funktion
 * @param b Eingabematrix
 * @param len Anzahl der Spalten bzw. Zeilen von b
 */
__global__ void kernel_computeB_sum(real_t *b, const unsigned int len) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len) {
    real_t sum = 0.0;
    for (unsigned int j = 0; j < len; ++j) {
      sum +=
          (j < tid) ? b[(tid * (tid + 1) / 2) + j] : b[(j * (j + 1) / 2) + tid];
    }
    b[(tid * (tid + 1) / 2) + tid] = -sum;
  }
}

extern "C" tm_t *cuda_computeBofZ(const tm_t *delta, const tm_t *dz,
                                  const tm_t *w) {
  tm_t *b = initTM(delta->size);
  real_t *dev_delta, *dev_dz, *dev_w = NULL, *dev_b;
  CUERROR(cudaMalloc((void **)&dev_b, b->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_delta, delta->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_dz, dz->num_elems * sizeof(real_t)));
  if (w != NULL) {
    CUERROR(cudaMalloc((void **)&dev_w, w->num_elems * sizeof(real_t)));
    CUERROR(cudaMemcpy(dev_w, w->elems, w->num_elems * sizeof(real_t),
                       cudaMemcpyHostToDevice));
  }
  CUERROR(cudaMemcpy(dev_delta, delta->elems, delta->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  CUERROR(cudaMemcpy(dev_dz, dz->elems, dz->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  const unsigned int gridsize_frac =
      (b->num_elems % THREADS_PER_BLOCK == 0)
          ? b->num_elems / THREADS_PER_BLOCK
          : (b->num_elems / THREADS_PER_BLOCK) + 1;
  kernel_computeB_frac<<<gridsize_frac, THREADS_PER_BLOCK>>>(
      dev_delta, dev_dz, dev_w, dev_b, b->num_elems, b->size);
  const unsigned int gridsize_sum = (b->size % THREADS_PER_BLOCK == 0)
                                        ? b->size / THREADS_PER_BLOCK
                                        : (b->size / THREADS_PER_BLOCK) + 1;
  kernel_computeB_sum<<<gridsize_sum, THREADS_PER_BLOCK>>>(dev_b, b->size);
  CUERROR(cudaDeviceSynchronize());
  CUERROR(cudaMemcpy(b->elems, dev_b, b->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_delta);
  cudaFree(dev_dz);
  cudaFree(dev_b);
  if (w != NULL)
    cudaFree(dev_w);
  return b;
}

extern "C" void cuda_computeBofZ_nomem(const real_t *dev_delta,
                                       const unsigned int delta_size,
                                       const unsigned int delta_num_elems,
                                       const real_t *dev_d, const real_t *dev_w,
                                       real_t *dev_bz) {
  const unsigned int gridsize_frac =
      (delta_num_elems % THREADS_PER_BLOCK == 0)
          ? delta_num_elems / THREADS_PER_BLOCK
          : (delta_num_elems / THREADS_PER_BLOCK) + 1;
  kernel_computeB_frac<<<gridsize_frac, THREADS_PER_BLOCK>>>(
      dev_delta, dev_d, dev_w, dev_bz, delta_num_elems, delta_size);
  const unsigned int gridsize_sum = (delta_size % THREADS_PER_BLOCK == 0)
                                        ? delta_size / THREADS_PER_BLOCK
                                        : (delta_size / THREADS_PER_BLOCK) + 1;
  kernel_computeB_sum<<<gridsize_sum, THREADS_PER_BLOCK>>>(dev_bz, delta_size);
}

__global__ void kernel_generateRandM(real_t *rand, int len) {}
