/*******************************************************************
 ***    poisson2D: Numerical solution of the Poisson PDE in 2D.
 ***
 ***    Solver function and iterators.
 ***
 *** Author: Nikos Tryfonidis, December 2015
 *** The MIT License (MIT)
 *** Copyright (c) 2015 Nikos Tryfonidis
 *** See LICENSE.txt
 *******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double residual(double **u, double **rhs, double dx, double dy, int nX, int nY);

/*  Poisson 2D Solver: Iterates until solution with desired accuracy is found 
   or until it reaches "maxIterations". Calls iterative functions from "iterative.c".
*/
double **solvePoisson2D(double ** (*iterate)(double **, double **, double, double, int, int), 
                        double **array2D, double **rhs, double dx, double dy, int nX, int nY)
{
    int t, nIterations, iterationsPerCheck, maxIterations;
    double res, tolerance;
    
    /* Set tolerance (value for residual considered adequate), 
       maximum iterations and iterations per residual check */
    tolerance = sqrt(nX*nY)*dx*dy; // Some multiple of expected accuracy of scheme
    maxIterations = 100*nX*nY; // Some multiple of expected number of iterations to converge
    iterationsPerCheck = 10;
    nIterations = 0; //counter

    /* Iterate until residual < tolerance or until maxIterations are reached */
    do {
        // Calculate residual
        res = residual(array2D, rhs, dx, dy, nX, nY);
        // Do "iterationsPerCheck" iterations before checking residual
        for (t=0;t<iterationsPerCheck;t++) {
            array2D = (*iterate)(array2D, rhs, dx, dy, nX, nY);
            nIterations += 1;
        }
        
        // Test print for residual every "iterationsPerCheck"
        //printf("iteration %d: residual = %f\ttolerance: %f\n", nIterations, res, tolerance);

    } while (res > tolerance && nIterations <= maxIterations);

    /* Print number of iterations needed */
    printf("Done! Iterations required: %d\n", nIterations);
    if(nIterations>=maxIterations) {
        printf("Warning: Maximum number of iterations reached! (%d)\n", nIterations);
    }

    return array2D;
}

/* Residual function: Calculates the residual
   after the application of relaxation.
*/
double residual(double **u, double **rhs, double dx, double dy, int nX, int nY) {
    int i, j;
    double Ax_r, res; 

    res = 0;
    /* Calculate "A*x - rho" and then add to res. */
    #pragma omp parallel for schedule(static) shared(u, rhs) private(i, j, Ax_r)\
                                firstprivate(nX, nY, dx, dy) reduction(+:res) default(none)
    for (i=1; i<nX-1; i++) {
        for (j=1; j<nY-1; j++) {    
            Ax_r = ( 1.0/(dx*dx)*(u[i-1][j]-2.0*u[i][j]+u[i+1][j]) +
                     1.0/(dy*dy)*(u[i][j-1]-2.0*u[i][j]+u[i][j+1]) -
                     rhs[i][j] );
            res += Ax_r*Ax_r;
        }
    }

    return sqrt(res);
}
