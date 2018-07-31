/******************************************************************************
 * Copyright (c) 2015 - 2016 Philipp Schubert.                                *
 * All rights reserved. This program and the accompanying materials are made  *
 * available under the terms of LICENSE.txt.                                  *
 *                                                                            *
 * Contributors:                                                              *
 *     Philipp Schubert                                                       *
 *****************************************************************************/

/** @file cuda_linalg.cu
 *  @brief Implementation der Prototypen aus cuda_linalg.cuh.
 *
 *  Hier sind die Wrapperfunktionen der Prototypen aus cuda_linalg.cuh
 *  implementiert, sowie entsprechende CUDA Kernel, die die Berechnungen auf der
 *  CUDA Karte durchführt.
 *
 *  @author Philipp D. Schubert
 *  @bug Keine Bugs bekannt.
 */

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>

// my c includes
extern "C" {
#include "m.h"
#include "tm.h"
#include "utils.cuh"
#include "utils.h"
}

/**
 * @brief Größe der TILES die von den Kernel Funktionen in diesem Modul benutzt
 * werden sollen.
 */
#define TILE_WIDTH 16
static const int THREADS_PER_BLOCK = 512;

/**
 * @brief Test Kernel zur Überprüfung des Moduls cuda_linalg.cu.
 */
__global__ void kernel_linalg_test() { printf("cuda_linalg: working fine!\n"); }

extern "C" void cuda_linalg_test() {
  kernel_linalg_test<<<1, 1>>>();
  cudaDeviceSynchronize();
}

/**
 * @brief Kernel Funktion zur Berechnung einer allgemeinen Matrixmultiplikation.
 *
 * Diese Kernelfunktion berechnet das Matrixprodukt zwischen zwei Matrizen von
 * beliebiger Dimension.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabematrix
 * @param arows Anzahl der Zeilen von a
 * @param acols Anzahl der Spalten von a
 * @param b Eingabematrix
 * @param brows Anzahl der Zeilen von b
 * @param bcols Anzahl der Spalten von b
 * @param splen Anzahl der Teilskalarprodukte, um das Gesamtergebnis einen
 * Skalarproduktes zu erhalten
 * @param c Ergebnismatrix
 */
__global__ void kernel_MmultM(const real_t *a, const unsigned int arows,
                              const unsigned int acols, const real_t *b,
                              const unsigned int brows,
                              const unsigned int bcols,
                              const unsigned int splen, real_t *c) {
  __shared__ real_t At[TILE_WIDTH][TILE_WIDTH];
  __shared__ real_t Bt[TILE_WIDTH][TILE_WIDTH];
  int bx = blockIdx.x;
  int by = blockIdx.y;
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int prow = bx * TILE_WIDTH + ty;
  int pcol = by * TILE_WIDTH + tx;
  real_t c_tmp = 0.0;
  for (int m = 0; m < splen; ++m) {
    /* int aax = prow; */
    int aay = m * TILE_WIDTH + tx;
    int bbx = m * TILE_WIDTH + ty;
    /* int bby = pcol; */
    At[ty][tx] = (prow < arows && aay < acols) ? a[prow * acols + aay] : 0.0;
    Bt[ty][tx] = (bbx < brows && pcol < bcols) ? b[bbx * bcols + pcol] : 0.0;
    __syncthreads();
    for (int k = 0; k < TILE_WIDTH; ++k)
      c_tmp += At[ty][k] * Bt[k][tx];
    __syncthreads();
    if (prow < arows && pcol < bcols)
      c[prow * bcols + pcol] = c_tmp;
  }
}

extern "C" m_t *cuda_MmultM(const m_t *A, const m_t *B) {
  CERROR(A->cols != B->rows, "bad dimension for matrix multiplication");

  m_t *C = initM(A->rows, B->cols);
  real_t *dev_a, *dev_b, *dev_c;
  CUERROR(cudaMalloc((void **)&dev_a, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_b, B->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_c, C->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_a, A->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  CUERROR(cudaMemcpy(dev_b, B->elems, B->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  int gdim1 = (A->rows % TILE_WIDTH == 0) ? A->rows / TILE_WIDTH
                                          : (A->rows / TILE_WIDTH) + 1;
  int gdim2 = (B->cols % TILE_WIDTH == 0) ? B->cols / TILE_WIDTH
                                          : (B->cols / TILE_WIDTH) + 1;
  int scalprodlen = (A->cols % TILE_WIDTH == 0) ? A->cols / TILE_WIDTH
                                                : (A->cols / TILE_WIDTH) + 1;
  dim3 gridsize(gdim1, gdim2);
  dim3 blocksize(TILE_WIDTH, TILE_WIDTH);
  kernel_MmultM<<<gridsize, blocksize>>>(dev_a, A->rows, A->cols, dev_b,
                                         B->rows, B->cols, scalprodlen, dev_c);
  CUERROR(cudaMemcpy(C->elems, dev_c, C->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
  return C;
}

/**
 * @brief Berechnet die Summe zweier Matrizen gleiche Dimension.
 *
 * Es wird die Summe zweier Matrizen gebildet, die die gleichen Dimensionen
 * besitzen.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabematrix
 * @param b Eingabematrix
 * @param len Anzahl der Matrixeinträge von a und b
 * @param c Ergebnismatrix
 */
__global__ void kernel_MaddM(const real_t *a, const real_t *b,
                             const unsigned int len, real_t *c) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len)
    c[tid] = a[tid] + b[tid];
}

extern "C" m_t *cuda_MaddM(const m_t *A, const m_t *B) {
  CERROR(A->rows != B->rows || A->cols != B->cols,
         "bad dimension for matrix multiplication");
  m_t *C = initM(A->rows, A->cols);
  real_t *dev_a, *dev_b, *dev_c;
  CUERROR(cudaMalloc((void **)&dev_a, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_b, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_c, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_a, A->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  CUERROR(cudaMemcpy(dev_b, B->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  int gridsize = (A->num_elems % THREADS_PER_BLOCK) == 0
                     ? A->num_elems / THREADS_PER_BLOCK
                     : (A->num_elems / THREADS_PER_BLOCK) + 1;
  kernel_MaddM<<<gridsize, THREADS_PER_BLOCK>>>(dev_a, dev_b, A->num_elems,
                                                dev_c);
  CUERROR(cudaMemcpy(C->elems, dev_c, A->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
  return C;
}

/**
 * @brief Berechnet das Vielfache einer beliebigen Matrix.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabematrix
 * @param x skalarer Wert
 * @param len Anzahl der Matrixeinträge von a
 * @param b Ergebnismatrix
 */
__global__ void kernel_Mmultscalar(const real_t *a, const real_t x,
                                   const unsigned int len, real_t *b) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len) {
    b[tid] = a[tid] * x;
  }
}

extern "C" m_t *cuda_Mmultscalar(const m_t *A, const real_t x) {
  m_t *B = initM(A->rows, A->cols);
  real_t *dev_a, *dev_b;
  CUERROR(cudaMalloc((void **)&dev_a, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_b, B->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_a, A->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  const unsigned int gridsize = (A->num_elems % THREADS_PER_BLOCK) == 0
                                    ? A->num_elems / THREADS_PER_BLOCK
                                    : (A->num_elems / THREADS_PER_BLOCK) + 1;
  kernel_Mmultscalar<<<gridsize, THREADS_PER_BLOCK>>>(dev_a, x, A->num_elems,
                                                      dev_b);
  CUERROR(cudaMemcpy(B->elems, dev_b, B->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_a);
  cudaFree(dev_b);
  return B;
}

/**
 * @brief Berechnet das Vielfache einer Matrix.
 *
 * Diese Funktion berechnet das Vielfache einer Matrix und schreibt das Ergebnis
 * direkt in die Eingabematrix zurück.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabe- / Ergebnismatrix
 * @param x skalarer Wert
 * @param len Anzahl der Einträge der Matrix a
 */
__global__ void kernel_Mmultscalar_noret(real_t *a, const real_t x,
                                         const unsigned int len) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len) {
    a[tid] = a[tid] * x;
  }
}

extern "C" void cuda_Mmultscalar_nomem(real_t *dev_a, const real_t x,
                                       const unsigned int len) {
  const unsigned int gridsize = (len % THREADS_PER_BLOCK) == 0
                                    ? len / THREADS_PER_BLOCK
                                    : (len / THREADS_PER_BLOCK) + 1;
  kernel_Mmultscalar_noret<<<gridsize, THREADS_PER_BLOCK>>>(dev_a, x, len);
}

/**
 * @brief Addiert einen skalaren Werte auf jeden Eintrag einer Matrix.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabematrix
 * @param x skalarer Wert
 * @param len Anzahl der Matrixeinträge
 * @param b Ergebnismatrix
 */
__global__ void kernel_Maddscalar(const real_t *a, const real_t x,
                                  const unsigned int len, real_t *b) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len)
    b[tid] = a[tid] + x;
}

extern "C" m_t *cuda_Maddscalar(const m_t *A, const real_t x) {
  m_t *B = initM(A->rows, A->cols);
  real_t *dev_a, *dev_b;
  CUERROR(cudaMalloc((void **)&dev_a, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_b, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_a, A->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  int gridsize = (A->num_elems % THREADS_PER_BLOCK) == 0
                     ? A->num_elems / THREADS_PER_BLOCK
                     : (A->num_elems / THREADS_PER_BLOCK) + 1;
  kernel_Maddscalar<<<gridsize, THREADS_PER_BLOCK>>>(dev_a, x, A->num_elems,
                                                     dev_b);
  CUERROR(cudaMemcpy(B->elems, dev_b, A->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_a);
  cudaFree(dev_b);
  return B;
}

/**
 * @brief Berechnet die Summe zweier Dreieckmatrizen gleicher Dimension.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabematrix
 * @param b Eingabematrix
 * @param len Anzahl der Matrixeinträge von a und b
 * @param c Ergebnismatrix
 */
__global__ void kernel_TMaddTM(const real_t *a, const real_t *b,
                               const unsigned int len, real_t *c) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len)
    c[tid] = a[tid] + b[tid];
}

extern "C" tm_t *cuda_TMaddTM(const tm_t *A, const tm_t *B) {
  CERROR(A->size != B->size, "bad dimension for matrix addition");
  tm_t *C = initTM(A->size);
  real_t *dev_a, *dev_b, *dev_c;
  CUERROR(cudaMalloc((void **)&dev_a, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_b, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_c, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_a, A->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  CUERROR(cudaMemcpy(dev_b, B->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  int gridsize = (A->num_elems % THREADS_PER_BLOCK) == 0
                     ? A->num_elems / THREADS_PER_BLOCK
                     : (A->num_elems / THREADS_PER_BLOCK) + 1;
  kernel_TMaddTM<<<gridsize, THREADS_PER_BLOCK>>>(dev_a, dev_b, A->num_elems,
                                                  dev_c);
  CUERROR(cudaMemcpy(C->elems, dev_c, A->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
  return C;
}

/**
 * @brief Berechnet das Vielfache einer Dreiecksmatrix.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabematrix
 * @param x skalarer Wert
 * @param len Anzahl der Matrixeinträge
 * @param b Ergebnismatrix
 */
__global__ void kernel_TMmultscalar(const real_t *a, const real_t x,
                                    const unsigned int len, real_t *b) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len)
    b[tid] = a[tid] * x;
}

extern "C" tm_t *cuda_TMmultscalar(const tm_t *A, const real_t x) {
  tm_t *B = initTM(A->size);
  real_t *dev_a, *dev_b;
  CUERROR(cudaMalloc((void **)&dev_a, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_b, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_a, A->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  int gridsize = (A->num_elems % THREADS_PER_BLOCK) == 0
                     ? A->num_elems / THREADS_PER_BLOCK
                     : (A->num_elems / THREADS_PER_BLOCK) + 1;
  kernel_TMmultscalar<<<gridsize, THREADS_PER_BLOCK>>>(dev_a, x, A->num_elems,
                                                       dev_b);
  CUERROR(cudaMemcpy(B->elems, dev_b, A->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_a);
  cudaFree(dev_b);
  return B;
}

/**
 * @brief Addiert einen skalaren Wert auf jeden Eintrag einer Matrix.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabematrix
 * @param x skalarer Wert
 * @param len Anzahl der Matrixeinträge
 * @param b Ergebnismatrix
 */
__global__ void kernel_TMaddscalar(const real_t *a, const real_t x,
                                   const unsigned int len, real_t *b) {
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len)
    b[tid] = a[tid] + x;
}

extern "C" tm_t *cuda_TMaddscalar(const tm_t *A, const real_t x) {
  tm_t *B = initTM(A->size);
  real_t *dev_a, *dev_b;
  CUERROR(cudaMalloc((void **)&dev_a, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_b, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_a, A->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  int gridsize = (A->num_elems % THREADS_PER_BLOCK) == 0
                     ? A->num_elems / THREADS_PER_BLOCK
                     : (A->num_elems / THREADS_PER_BLOCK) + 1;
  kernel_TMaddscalar<<<gridsize, THREADS_PER_BLOCK>>>(dev_a, x, A->num_elems,
                                                      dev_b);
  CUERROR(cudaMemcpy(B->elems, dev_b, A->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_a);
  cudaFree(dev_b);
  return B;
}

/**
 * @brief Berechnet das Matrixprodukt aus einer Dreiecks- und einer allgemeinen
 * Matrix.
 *
 * Es wird das Matrixprodukt aus einer Dreiecks- und einer allgemeinen Matrix
 * gebildet. Es wird geteilter Speicher genutzt, um mehrfach benutzte Daten zu
 * cachen. Die Dreiecksmatrix wird dazu in den geteilten Speicher At geladen,
 * dazu muss ein bedingtes Lesen der Einträge aus Matrix a durchgeführt werden!
 * Die Einträge aus der Matrix b können ohne weiteres in den geteilten Speicher
 * Bt geladen werden.
 *
 * @warning GPU Kernel Funktion
 * @param a Eingabematrix - Dreiecksmatrix
 * @param asize Dimensionen der Matrix a
 * @param b Eingabematrix - allgemeine Matrix
 * @param brows Anzahl der Zeilen von b
 * @param bcols Anzahl der Spalten von b
 * @param splen Anzahl der Blöcke die benötigt werden, um ein Skalarprodukt
 * zwischen a und b zu bilden
 * @param c Ergebnismatrix
 */
__global__ void kernel_TMmultM(const real_t *a, const unsigned int asize,
                               const real_t *b, const unsigned int brows,
                               const unsigned int bcols,
                               const unsigned int splen, real_t *c) {
  __shared__ real_t At[TILE_WIDTH][TILE_WIDTH];
  __shared__ real_t Bt[TILE_WIDTH][TILE_WIDTH];
  int bx = blockIdx.x;
  int by = blockIdx.y;
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int prow = bx * TILE_WIDTH + ty;
  int pcol = by * TILE_WIDTH + tx;
  real_t c_tmp = 0.0;
  for (int m = 0; m < splen; ++m) {
    /* int aax = prow; */
    int aay = m * TILE_WIDTH + tx;
    int bbx = m * TILE_WIDTH + ty;
    /* int bby = pcol; */
    if (aay <= prow) {
      At[ty][tx] = (prow < asize && aay < asize)
                       ? a[(prow * (prow + 1)) / 2 + aay]
                       : 0.0;
    } else {
      At[ty][tx] =
          (prow < asize && aay < asize) ? a[(aay * (aay + 1)) / 2 + prow] : 0.0;
    }
    Bt[ty][tx] = (bbx < brows && pcol < bcols) ? b[bbx * bcols + pcol] : 0.0;
    __syncthreads();
    for (int k = 0; k < TILE_WIDTH; ++k)
      c_tmp += At[ty][k] * Bt[k][tx];
    __syncthreads();
    if (prow < asize && pcol < bcols)
      c[prow * bcols + pcol] = c_tmp;
  }
}

extern "C" m_t *cuda_TMmultM(const tm_t *A, const m_t *B) {
  CERROR(A->size != B->rows, "bad dimension for matrix multiplication");

  m_t *C = initM(A->size, B->cols);
  real_t *dev_a, *dev_b, *dev_c;
  CUERROR(cudaMalloc((void **)&dev_a, A->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_b, B->num_elems * sizeof(real_t)));
  CUERROR(cudaMalloc((void **)&dev_c, C->num_elems * sizeof(real_t)));
  CUERROR(cudaMemcpy(dev_a, A->elems, A->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  CUERROR(cudaMemcpy(dev_b, B->elems, B->num_elems * sizeof(real_t),
                     cudaMemcpyHostToDevice));
  const unsigned int gdim1 = (A->size % TILE_WIDTH == 0)
                                 ? A->size / TILE_WIDTH
                                 : (A->size / TILE_WIDTH) + 1;
  const unsigned int gdim2 = (B->cols % TILE_WIDTH == 0)
                                 ? B->cols / TILE_WIDTH
                                 : (B->cols / TILE_WIDTH) + 1;
  dim3 gridsize(gdim1, gdim2);
  dim3 blocksize(TILE_WIDTH, TILE_WIDTH);
  kernel_TMmultM<<<gridsize, blocksize>>>(dev_a, A->size, dev_b, B->rows,
                                          B->cols, gdim1, dev_c);
  CUERROR(cudaMemcpy(C->elems, dev_c, C->num_elems * sizeof(real_t),
                     cudaMemcpyDeviceToHost));
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
  return C;
}

extern "C" void cuda_TMmultM_nomem(const real_t *dev_a,
                                   const unsigned int asize,
                                   const real_t *dev_b,
                                   const unsigned int brows,
                                   const unsigned int bcols, real_t *dev_c) {
  CERROR(asize != brows, "bad dimension for matrix multiplication");
  const unsigned int gdim1 =
      (asize % TILE_WIDTH == 0) ? asize / TILE_WIDTH : (asize / TILE_WIDTH) + 1;
  const unsigned int gdim2 =
      (bcols % TILE_WIDTH == 0) ? bcols / TILE_WIDTH : (bcols / TILE_WIDTH) + 1;
  dim3 gridsize(gdim1, gdim2);
  dim3 blocksize(TILE_WIDTH, TILE_WIDTH);
  kernel_TMmultM<<<gridsize, blocksize>>>(dev_a, asize, dev_b, brows, bcols,
                                          gdim1, dev_c);
}
