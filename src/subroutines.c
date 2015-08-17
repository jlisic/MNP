
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R.h> 
#include <R_ext/Lapack.h>
#include "vector.h"
#include "rand.h"

/*  The Sweep operator */
/*
 * Goodnight, American Statistician 1979 is a good reference.
 * If you also happen to have Little and Rubin's book "statistical 
 * analysis with missing data" handy, it's in there on page 148.
 *
 * example for the 3x3 matrix A:
 *
 * a11 a12 a13
 * a21 a22 a23
 * a31 a32 a33
 *
 * Sweep(A,1)
 *
 * -1/a11     | a12/a11            a13/a11
 *  ----------|--------------------------------------
 *            | a22 - a12*a12/a11  a23 - a13*a12/a11 
 *            |                    a33 - a13*a13/a11 
 *
 * Note: due to symmetry (we are only working with square symmetric matricies here
 *       only the upper triangular matrix is calucluated
 *
 * Note: there are a few different versions floating around in the literature.
 * Note (jjl): IIRC there are some things that can be done for symmetric matricies
 *   but I forgot what they where :( 
 */
void SWP(
	 double **X,             /* The Matrix to work on */
	 int k,                  /* The row to sweep */
	 int size)               /* The dim. of X */
{
  int i,j;

  if (X[k][k] < 10e-20) 
    error("SWP: singular matrix.\n");
  else
    X[k][k]=-1/X[k][k];
  for(i=0;i<size;i++)
    if(i!=k){
      X[i][k]=-X[i][k]*X[k][k];
      X[k][i]=X[i][k];
    }
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      if(i!=k && j!=k)
	X[i][j]=X[i][j]+X[i][k]*X[k][j]/X[k][k];
  
}

/* R interface */
void RSWP( 
   double * x,
   int * kPtr,
   int * sizePtr
   ) {

  size_t i; 
  double ** X;

  X = calloc( *sizePtr, sizeof(double *) );
  for(i = 0; i < *sizePtr; i++) X[i] = &x[ *sizePtr * i ];

  SWP(X, *kPtr, *sizePtr);

  free(X);
  return;
} 



/* inverting a matrix */
void dinv(
  double **X,
  int	size,
  double **X_inv
  ) {
  int i,j, k, errorM;
  double *pdInv = doubleArray(size*(size+1)/2);

  for (i = 0, j = 0; j < size; j++) 
    for (k = 0; k <= j; k++) pdInv[i++] = X[k][j];

  /* ref: http://www.netlib.org/lapack/explore-html/d2/d45/dpptrf_8f.html */
  /* 
   * dpprtf:
   * Ax = b
   * U     - Upper triangle of A is stored
   * &size - order of the matrix A 
   * pdInv - in/out, using N*N_1)/2 elements (N is an int).
   *   e.g.
   *         a11 a12 a13 a14
   *             a22 a23 a24
   *                 a33 a34
   *                     a44
   *
   * error = 0 successful
   *       < 0 (-i) the "ith" value had an illegal value
   *       > 0 (i) the leading minor of order i is not PD
   */
  F77_CALL(dpptrf)("U", &size, pdInv, &errorM);  /* cholesky factorization via lapack */
  if (!errorM) {
    F77_CALL(dpptri)("U", &size, pdInv, &errorM);
    if (errorM) {
      Rprintf("LAPACK dpptri failed, %d\n", errorM);
      error("Exiting from dinv().\n");
    }
  }
  else {
    Rprintf("LAPACK dpptrf failed, %d\n", errorM);
    error("Exiting from dinv().\n");
  }
  /* I like this code :)
   *
   * 0 1 2 3 
   *   5 6 7
   *     8 9
   *       10
   */
  for (i = 0, j = 0; j < size; j++) {
    for (k = 0; k <= j; k++) {
      X_inv[j][k] = pdInv[i];
      X_inv[k][j] = pdInv[i++];
    }
  }

  free(pdInv);
}


/* Cholesky decomposition */
/* returns lower triangular matrix */
void dcholdc(double **X, int size, double **L)
{
  int i, j, k, errorM;
  double *pdTemp = doubleArray(size*size);

  for (j = 0, i = 0; j < size; j++) 
    for (k = 0; k <= j; k++) 
      pdTemp[i++] = X[k][j];
  /* ref: http://www.netlib.org/lapack/explore-html/d2/d45/dpptrf_8f.html */
  /* 
   * dpprtf:
   * Ax = b
   * U     - Upper triangle of A is stored
   * &size - order of the matrix A 
   * pdInv - in/out, using N*N_1)/2 elements (N is an int).
   *   e.g.
   *         a11 a12 a13 a14
   *             a22 a23 a24
   *                 a33 a34
   *                     a44
   *
   * error = 0 successful
   *       < 0 (-i) the "ith" value had an illegal value
   *       > 0 (i) the leading minor of order i is not PD
   */
  F77_CALL(dpptrf)("U", &size, pdTemp, &errorM);
  if (errorM) {
    Rprintf("LAPACK dpptrf failed, %d\n", errorM);
    error("Exiting from dcholdc().\n");
  }
  for (j = 0, i = 0; j < size; j++) {
    for (k = 0; k < size; k++) {
      if(j<k)
	L[j][k] = 0.0;
      else
	L[j][k] = pdTemp[i++];
    }
  }

  free(pdTemp);
} 

/* calculate the determinant of the positive definite symmetric matrix
   using the Cholesky decomposition  */
double ddet(double **X, int size, int give_log)
{
  int i;
  double logdet=0.0;
  double **pdTemp = doubleMatrix(size, size);
  
  dcholdc(X, size, pdTemp);
  for(i = 0; i < size; i++)
    logdet += log(pdTemp[i][i]);

  FreeMatrix(pdTemp, size);
  if(give_log)
    return(2.0*logdet);
  else
    return(exp(2.0*logdet));
}
