/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/** @file cuda_disfuncs.cu
 *  @brief Implementation der Prototypen aus cuda_disfuncs.cuh.
 *
 *  Hier sind die Wrapperfunktionen der Prototypen aus cuda_disfuncs.cuh
 *  implementiert, sowie entsprechende CUDA Kernel, die die Berechnungen auf der
 *  CUDA Karte durchführt.
 *
 *  @author Philipp D. Schubert
 *  @bug Keine Bugs bekannt.
 */

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Zu benutzende Größe der TILES im den Kernels dieses Moduls.
 */
#define TILE_WIDTH 16

// my C includes
extern "C" {
#include "m.h"
#include "tm.h"
#include "utils.cuh"
#include "utils.h"
}

/**
 * @brief Test Kernel zur Überprüfung des Moduls cuda_disfuncs.cu.
 */
__global__ void kernel_disfuncs_test() {
  printf("cuda_disfuncs: working fine!\n");
}

extern "C" void cuda_disfuncs_test() {
  kernel_disfuncs_test<<<1, 1>>>();
  cudaDeviceSynchronize();
}

/**
 * @brief Berechnet aus gegebener Datenmatrix eine Distanzmatrix der
 * euklidischen Distanzen.
 *
 * Dieser Kernel berechnet aus gegebener Datenmatrix einer Matrix, die die
 * euklidischen Distanzen enthält. Zur Berechnung dieser Distanzmatrix wird eine
 * spezielle Variante der allgemeinen Matrixmultiplikation aus cuda_linalg.cu
 * verwendet. Der geteilte Speicher At und Bt wird dazu ausgenutzt. At wird mit
 * Daten aus der Matrix a befüllt, Bt mit den transponierten Daten aus der
 * Matrix a. Auf diese weise kann die Distanzmatrix sehr geschickt berechnet
 * werden, wobei die Graphikkarte nahezu optimal ausgelastet wird.
 *
 * @warning GPU Kernel Funktion
 * @param a Datenmatrix
 * @param arows Anzahl der Zeilen der Datenmatrix
 * @param acols Anzahl der Spalten der Datenmatrix
 * @param splen Anzahl der mindestens benötigten Blöcke, um die Matrix in
 * Dimension 1 abzudecken
 * @param dist untere Dreiecksmatrix in die die Distanzen geschrieben werden
 */
__global__ void kernel_distanceMatrix(const real_t *a, const unsigned int arows,
                                      const unsigned int acols,
                                      const unsigned int splen, real_t *dist) {
  __shared__ real_t At[TILE_WIDTH][TILE_WIDTH];
  __shared__ real_t Bt[TILE_WIDTH][TILE_WIDTH];
  int bx = blockIdx.x;
  int by = blockIdx.y;
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int prow = bx * TILE_WIDTH + ty;
  int pcol = by * TILE_WIDTH + tx;
  real_t dist_tmp = 0.0;
  real_t diff;
  for (int m = 0; m < splen; ++m) {
    /* int aax = prow; */
    int aay = m * TILE_WIDTH + tx;
    int bbx = m * TILE_WIDTH + ty;
    /* int bby = pcol; */
    At[ty][tx] = (prow < arows && aay < acols) ? a[prow * acols + aay] : 0.0;
    Bt[ty][tx] = (bbx < acols && pcol < arows) ? a[pcol * acols + bbx] : 0.0;
    __syncthreads();
    for (int k = 0; k < TILE_WIDTH; ++k) {
      diff = (At[ty][k] - Bt[k][tx]);
      dist_tmp += diff * diff;
    }
    __syncthreads();
    if (prow < arows && pcol < arows && (pcol <= prow))
      dist[(prow * (prow + 1)) / 2 + pcol] = sqrtf(dist_tmp);
  }
}

extern "C" tm_t *cuda_distanceMatrix(const m_t *data) {
  tm_t *dist = initTM(data->rows);
  real_t *dev_data, *dev_dist;
  CUERROR(cudaMalloc((void **)&dev_data, data->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_dist, dist->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_data, data->elems, data->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  int gdim1, gdim2;
  gdim1 = gdim2 = (data->rows % TILE_WIDTH == 0)
                      ? data->rows / TILE_WIDTH
                      : (data->rows / TILE_WIDTH) + 1;
  int scalprodlen = (data->cols % TILE_WIDTH == 0)
                        ? data->cols / TILE_WIDTH
                        : (data->cols / TILE_WIDTH) + 1;
  dim3 gridsize(gdim1, gdim2);
  dim3 blocksize(TILE_WIDTH, TILE_WIDTH);
  kernel_distanceMatrix<<<gridsize, blocksize>>>(
      dev_data, data->rows, data->cols, scalprodlen, dev_dist);
  CUERROR(cudaMemcpy(dist->elems, dev_dist, dist->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_data);
  cudaFree(dev_dist);
  return dist;
}

extern "C" void cuda_distanceMatrix_nomem(const real_t *dev_x,
                                          unsigned int xrows,
                                          unsigned int xcols, real_t *dev_d) {
  int gdim1, gdim2;
  gdim1 = gdim2 =
      (xrows % TILE_WIDTH == 0) ? xrows / TILE_WIDTH : (xrows / TILE_WIDTH) + 1;
  int scalprodlen =
      (xcols % TILE_WIDTH == 0) ? xcols / TILE_WIDTH : (xcols / TILE_WIDTH) + 1;
  dim3 gridsize(gdim1, gdim2);
  dim3 blocksize(TILE_WIDTH, TILE_WIDTH);
  kernel_distanceMatrix<<<gridsize, blocksize>>>(dev_x, xrows, xcols,
                                                 scalprodlen, dev_d);
}
