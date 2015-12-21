/*******************************************************************
 ***    poisson2D: Numerical solution of the Poisson PDE in 2D
 ***
 ***    Solver function and iterators.
 ***
 ***    Author: Nikos Tryfonidis
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../headers/memory.h"
#include "../headers/definitions.h"

/*******************************************
 Iterative Methods here:
 *******************************************/

/* - Jacobi
   Performs a single Jacobi iteration over the 2D domain.
   Does NOT assume that dx = dy.
 */
double **jacobiIteration2D(double **u_new2D, double **u, double **f, 
                            int nX, int nY, double dx, double dy)
{
    int i, j;

    for(i=1;i<nX-1;i++) {
        for(j=1;j<nY-1;j++) {

            u_new2D[i][j] = ((dy*dy)*(u[i-1][j] + u[i+1][j]) + 
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy)) );
        }
    }

    return u_new2D;
}

/* - Gauss-Seidel
   Performs a single Gauss-Seidel iteration over the 2D domain.
   Does NOT assume that dx = dy.
 */
double **gaussSeidelIteration2D(double **u, double **f, 
                            int nX, int nY, double dx, double dy)
{
    int i, j;

    for(i=1;i<nX-1;i++) {
        for(j=1;j<nY-1;j++) {

            u[i][j] = ((dy*dy)*(u[i-1][j] + u[i+1][j]) + 
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy)) );
        }
    }

    return u;
}


/********************************************
 Iterative Methods End
 ********************************************/

/* Copies array and returns */
double **copyArray2D(double **a, double **from, int nX, int nY) {
    int i, j;

    for(i=0;i<nX;i++) {
        for(j=0;j<nY;j++) {
            a[i][j] = from[i][j];
        }
    }

    return a;
}

/* Residual function: Calculates the residual
   after the application of relaxation
*/
double residual(double **u, double **rhs, double dx, double dy, int nX, int nY) {
    int i, j;
    double Ax_r, res; 

    res = 0;
    /* Calculate "A*x - rho" and then add to res. */
    for (i=1; i<nX-1; i++) {
        for (j=1; j<nY-1; j++) {    
            Ax_r = ( 1.0/(dx*dx)*(u[i-1][j]-2.0*u[i][j]+u[i+1][j]) +
                     1.0/(dy*dy)*(u[i][j-1]-2.0*u[i][j]+u[i][j+1]) -
                     rhs[i][j] );
            res += sqrt(Ax_r*Ax_r);
        }
    }

    return res;
}

/* Iterates until solution with desired accuracy is found */
double **solvePoisson2D(double **array2D, double **rhs, double dx, double dy, int nX, int nY)
{
    int i, j, t, nIterations, iterationsPerCheck;
    double res, tolerance;
    double **array2D_old;
    
    // Allocate 2D array to store previous timestep solution
    array2D_old = array2D_contiguous (nX, nY);

    /* Iterate until residual < tolerance */
    tolerance = nX*nY*pow(10, -6);
    nIterations = 0;
    iterationsPerCheck = 10;
    do {
        // Do "iterationsPerCheck" iterations before checking residual
        for (t=0;t<iterationsPerCheck;t++) {
            // Copy to "old" (needed for Jacobi)
            //copyArray2D(array2D_old, array2D, nX, nY);

            // Do one Jacobi iteration 
            //array2D = jacobiIteration2D(array2D, array2D_old, rhs, nX, nY, dx, dy);

            // Do one Gauss-Seidel iteration 
            array2D = gaussSeidelIteration2D(array2D, rhs, nX, nY, dx, dy);

            nIterations += 1;
        }
        // Calculate residual
        res = residual(array2D, rhs, dx, dy, nX, nY);

        // test-remove
        printf("iteration %d: residual = %f\ttolerance: %f\n", nIterations, res, tolerance);

        
    } while (res > tolerance);

    printf("Iterations: %d\n", nIterations);

    free_array2D_contiguous(array2D_old);
    return array2D;
}
