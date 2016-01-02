/*******************************************************************
 ***    poisson2D: Numerical solution of the Poisson PDE in 2D.
 ***    Iterative method functions, called by solvePoisson2D function
 ***    in poisson2D.c
 ***
 *** Author: Nikos Tryfonidis, December 2015
 *** The MIT License (MIT)
 *** Copyright (c) 2015 Nikos Tryfonidis
 *** See LICENSE.txt
 *******************************************************************/

#include <math.h>
#include <omp.h>

#include "../headers/memory.h"

double **copyArray2D(double **a, double **from, int nX, int nY); //Needed for Jacobi 


/* - Jacobi Method.
   Performs a single Jacobi iteration over the 2D domain.

   The decision to include "old" u allocation and copy inside
   the iterator function was made, in order to keep the same form
   with the other iterators.
   Performance was not found to suffer significantly from this.

   Does NOT assume that dx = dy.
 */
double **jacobiIteration2D(double **u, double **f, 
                           double dx, double dy, int nX, int nY)
{
    int i, j;
    double **u_old;
    u_old = array2D_contiguous (nX, nY);
    u_old = copyArray2D(u_old, u, nX, nY); 

    for(i=1;i<nX-1;i++) {
        for(j=1;j<nY-1;j++) {

            u[i][j] = ((dy*dy)*(u_old[i-1][j] + u_old[i+1][j]) + 
                             (dx*dx)*(u_old[i][j-1] + u_old[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy)) );
        }
    }

    free_array2D_contiguous(u_old);
    return u;
}

/* - Weighted Jacobi Method.
   Performs a single weighted Jacobi iteration over the 2D domain.

   Other than the given w factor, it is similar to the simple 
   Jacoi method.

   Does NOT assume that dx = dy.
 */
double **weightedJacobiIteration2D(double **u, double **f, 
                           double dx, double dy, int nX, int nY)
{
    int i, j;
    double w = 0.5;
    double **u_old, u_tmp;
    u_old = array2D_contiguous (nX, nY);
    u_old = copyArray2D(u_old, u, nX, nY); 

    for(i=1;i<nX-1;i++) {
        for(j=1;j<nY-1;j++) {
            u_tmp = ((dy*dy)*(u_old[i-1][j] + u_old[i+1][j]) + 
                             (dx*dx)*(u_old[i][j-1] + u_old[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy)) );
            u[i][j] = u_old[i][j] + w*(u_tmp - u_old[i][j]);
        }
    }

    free_array2D_contiguous(u_old);
    return u;
}


/* - Gauss-Seidel
   Performs a single Gauss-Seidel iteration over the 2D domain.
   Does NOT assume that dx = dy.
 */
double **gaussSeidelIteration2D(double **u, double **f, 
                           double dx, double dy, int nX, int nY)
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

/* - Red-Black Gauss-Seidel
   Performs a single Red-Black Gauss-Seidel iteration over the 2D domain.
   Does NOT assume that dx = dy.
 */
double **redBlackGaussSeidelIteration2D(double **u, double **f, 
                                double dx, double dy, int nX, int nY)
{
    int i, j;

    /* RED sweep (odd -> odd, even -> even ) */
    for(i=1;i<nX-1;i++) {
        for(j=(i%2)+1;j<nY-1;j+=2) {

            u[i][j] = ((dy*dy)*(u[i-1][j] + u[i+1][j]) + 
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy)) );
        }
    }

    /* BLACK sweep (odd -> even, even -> odd ) */
    for(i=1;i<nX-1;i++) {
        for(j=((i+1)%2) + 1;j<nY-1;j+=2) {

            u[i][j] = ((dy*dy)*(u[i-1][j] + u[i+1][j]) + 
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy)) );
        }
    }

    return u;
}

/* - Successive Overrelaxation (SOR)
   Performs a single SOR iteration over the 2D domain.
   Value for w (overrelaxation factor) is set inside the function.
   From R. Leveque, ch. 4.2.2: optimal w for Poisson equation
   is 2-2*pi*h.

   Does NOT assume that dx = dy.
 */
double **SORIteration2D(double **u, double **f, 
                           double dx, double dy, int nX, int nY)
{
    int i, j;
    
    double w = 2.0-2.0*M_PI*dx;

    for(i=1;i<nX-1;i++) {
        for(j=1;j<nY-1;j++) {

            u[i][j] = w*((dy*dy)*(u[i-1][j] + u[i+1][j]) + 
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy))) +
                      (1.0-w)*u[i][j];
        }
    }

    return u;
}

/* - Red - Black Successive Overrelaxation (SOR)
   Performs a single SOR iteration over the 2D domain.
   Value for w (overrelaxation factor) is set inside the function.
   From R. Leveque, ch. 4.2.2: optimal w for Poisson equation
   is 2-2*pi*h.

   Does NOT assume that dx = dy.
 */
double **redBlackSORIteration2D(double **u, double **f, 
                           double dx, double dy, int nX, int nY)
{
    int i, j;
    
    double w = 2.0-2.0*M_PI*dx;

    /* RED sweep (odd -> odd, even -> even ) */
    for(i=1;i<nX-1;i++) {
        for(j=(i%2)+1;j<nY-1;j+=2) {

            u[i][j] = w*((dy*dy)*(u[i-1][j] + u[i+1][j]) + 
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy))) +
                      (1.0-w)*u[i][j];
        }
    }

    /* BLACK sweep (odd -> even, even -> odd ) */
    for(i=1;i<nX-1;i++) {
        for(j=((i+1)%2) + 1;j<nY-1;j+=2) {

            u[i][j] = w*((dy*dy)*(u[i-1][j] + u[i+1][j]) + 
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) - 
                             (dx*dx)*(dy*dy)*f[i][j])/(2.0*((dx*dx) + (dy*dy))) +
                      (1.0-w)*u[i][j];
        }
    }

    return u;
}

/* Copies array and returns (needed for Jacobi) */
double **copyArray2D(double **a, double **from, int nX, int nY) {
    int i, j;
    for(i=0;i<nX;i++) {
        for(j=0;j<nY;j++) {
            a[i][j] = from[i][j];
        }
    }

    return a;
}
