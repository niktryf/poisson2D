/*******************************************************************
 ***    poisson2D: Numerical solution of the Poisson PDE in 2D.
 ***    Solves (d^2/dx^2)u + (d^2/dy^2)u = f(x,y)
 ***    with given boundary conditions.
 ***
 ***    See README.txt for details
 ***
 *** Author: Nikos Tryfonidis, December 2015
 *** The MIT License (MIT)
 *** Copyright (c) 2015 Nikos Tryfonidis
 *** See LICENSE.txt
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "../headers/memory.h"
#include "../headers/setup.h"
#include "../headers/poisson2D.h"
#include "../headers/iterative.h"
#include "../headers/io.h"

/* Definitions on space (x,y range) */
#include "../headers/definitions.h"

int main(int argc, char *argv[]) {
    int nX, nY, i, j;
    double **u, **rhs;
    double **(*iterativeMethod)(double **, double **, double, double, int, int);
    double dx, dy;
    double t1, t2;

    /* Get command line arguments (nX, nY) */
    getArgs(argc, argv, &nX, &nY);

    // Calculate dx, dy (X0, XL, Y0, YL are defined in "definitions.h")
    dx = (XL - X0)/(nX-1);
    dy = (YL - Y0)/(nY-1);

    /* Allocate 2D arrays for u, rhs and initialize to zero */
    u = array2D_contiguous(nX, nY);
    rhs = array2D_contiguous(nX, nY);

    /* Set Boundary Conditions */
    u = setBoundaries2D(u, dx, dy, nX, nY);

    /* RHS */
    rhs = setRHS2D(rhs, dx, dy, nX, nY);

    /**************************************************
     *** Choose iterative method from "iterative.h" ***
     **************************************************/
    iterativeMethod = &gaussSeidelIteration2D;

    /**************************************************/

    /* Call Solver (and measure execution time) */
    printf("Iterating...  \n");
    t1 = omp_get_wtime();

    u = solvePoisson2D(iterativeMethod, u, rhs, dx, dy, nX, nY);

    t2 = omp_get_wtime();

    /* Write Grid Data and Output to files */
    writeGridData(X0, XL, Y0, YL, dx, dy, nX, nY, "output/gridData.txt");
    writeFileOutput(u, nX, nY, "output/output.txt");

    /* Report execution time */
    printf("Execution Time: %f seconds.\n", t2-t1);

    /* Free Memory */
    free_array2D_contiguous(u);
    free_array2D_contiguous(rhs);

    return 0;
}
