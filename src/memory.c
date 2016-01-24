/*******************************************************************
 ***    poisson2D: Numerical solution of the Poisson PDE in 2D.
 ***    Memory allocation functions for dynamic 2D arrays,
 ***    contiguous in memory.
 ***
 *** Author: Nikos Tryfonidis, December 2015
 *** The MIT License (MIT)
 *** Copyright (c) 2015 Nikos Tryfonidis
 *** See LICENSE.txt
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/* Allocates 2D dynamic array of doubles 
   with given dimensions, contiguous in memory. 
   Array is initialized to zeros.

   Note: initialization is done in parallel with OpenMP,
   mainly to achieve better cache utilization (and NUMA affinity)
   when calling this function in the Jacobi method.
*/
double **array2D_contiguous (int nRows, int nCols) 
{
  int i, j;
  double **array2D;

  //Allocate row storage:
  array2D = (double **)malloc(nRows * sizeof(double *));
  if (array2D == NULL) {
    printf("Out of memory (array2D_contiguous, pointer array)! Exiting...\n");
    exit(-1);
  }

  //Allocate big chunk o' memory:
  array2D[0] = (double *)malloc(nRows * nCols * sizeof(double));
  if (array2D[0] == NULL) {
    printf("Out of memory (array2D_contiguous, chunk array)! Exiting...\n");
    exit(-1);
  }

  //Point
  for (i=1;i<nRows;i++) {
    array2D[i] = array2D[i-1] + nCols;
  }

  // Initialize to zero:
  #pragma omp parallel for schedule(static) shared(array2D) private(i, j)\
                           firstprivate(nRows, nCols) default(none)
  for (i=0;i<nRows;i++) {
    for (j=0;j<nCols;j++) {
      array2D[i][j] = 0;
    }
  }

  return array2D;
}

/* Frees memory reserved for 2D array that has been created with
   the "array2D_contiguous" function.
 */
void free_array2D_contiguous(double **array2D)
{
  free(array2D[0]); 
  free(array2D);
}
