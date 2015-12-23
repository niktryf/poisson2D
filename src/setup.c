/*******************************************************************
 ***    poisson2D: Numerical solution of the Poisson PDE in 2D
 ***
 ***    Setup functions for Boundary Conditions and
 ***    Right-Hand-Side (RHS) function. 
 ***
 ***    Author: Nikos Tryfonidis
 *******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../headers/definitions.h"

/*******************************************************************************
 *** Boundary Conditions Function
 *******************************************************************************/

/* Function for Boundary Conditions fBC(x, y).
   Set desired function here:
*/
double fBC(double x, double y) 
{
    return sin(2.0*M_PI*(x+y));
}

/*  Boundary Conditions setup: Set up your bounds here!
    Sets bounds in 2D: (top-bottom, left-right)
    
                 t t t
                l     r
                l     r
                l     r
                 b b b  

    This function uses the fBC function to set the desired 
    Boundary Conditions.
*/
double **setBoundaries2D(double **array2D, double dx, double dy, int nX, int nY) 
{
    int i, j;
    /* For cases where BCs are a function of x, y: */
    double x, y;

    /* Top Boundary (row 0, all columns) */
    for(j=0;j<nY;j++) {
        x = 0;
        y = j*dy;
        array2D[0][j] = fBC(x, y);
    }

    /* Bottom Boundary (row nX-1, all columns) */
    for(j=0;j<nY;j++) {
        x = (nX-1)*dx;
        y = j*dy;
        array2D[nX-1][j] = fBC(x, y);
    }

    /* Left Boundary (column 0, all rows) */
    for(i=0;i<nX;i++) {
        x = i*dx;
        y = 0;
        array2D[i][0] = fBC(x, y);
    }

    /* Right Boundary (column nY-1, all rows) */
    for(i=0;i<nX;i++) {
        x = i*dx;
        y = (nY-1)*dy;
        array2D[i][nY-1] = fBC(x, y);
    }

    return array2D;
}

/*******************************************************************************
 *** Source (Right-Hand-Side) Function
 *******************************************************************************/

/* Right-Hand-Side function, f(x,y) */
double f(double x, double y) 
{
    return -8.0*M_PI*M_PI*sin(2.0*M_PI*(x+y));
}

/* Set the Right-Hand-Side (source term) of the Poisson PDE */
double **setRHS2D(double **array2D, double dx, double dy, int nX, int nY) 
{
    int i, j;

    for(i=1;i<nX-1;i++) {
        for(j=1;j<nY-1;j++) {
            array2D[i][j] = f(i*dx, j*dy);
        }
    }

    return array2D;
}
