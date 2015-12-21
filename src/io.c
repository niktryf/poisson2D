/*******************************************************************
 ***    poisson2D: Numerical solution of the Poisson PDE in 2D
 ***
 ***    Input and Output functions.
 ***
 ***    Author: Nikos Tryfonidis
 *******************************************************************/
#include <stdlib.h>
#include <stdio.h>

/*
  Reads arguments and sets timestep, length and output interval from command line
*/
int getArgs(int argc, char **argv, int *nX, int *nY)
{
  if (argc == 3) {
    sscanf(argv[1],"%d",&(*nX) );
    sscanf(argv[2],"%d",&(*nY) );

  }
  else {
    printf("\n*** User Error: Command line arguments needed! ***\n\n");
    printf("Please run with the following arguments:\n\n");
    printf("./poisson2D <number of data points in X> <number of data points in Y>\n\n");
    printf("Exiting...\n\n");
    exit(-1);
  }
  return 0;
}

/* Writes Grid data to file 
   Useful for analyzing results, visualization etc
*/
void writeGridData(double X0, double XL, double Y0, double YL, 
                    double dx, double dy, 
                    int nX, int nY, char *filename_GridData)
{
    FILE *outputFile;

    /* Open Output File */
    outputFile = fopen(filename_GridData, "w");

    /* Write Data */
    fprintf(outputFile, "# X - Y domain: Xstart Xend Ystart Yend : \n%f %f %f %f\n", X0, XL, Y0, YL);
    fprintf(outputFile, "# step size dx, dy: \n%f %f\n", dx, dy);
    fprintf(outputFile, "# number of grid points in X and Y dimensions: \n%d %d\n", nX, nY);

    /* Close Output File */
    fclose(outputFile);
}

/* Writes 2D solution array to file */
void writeFileOutput (double **outputArray, int nX, int nY, char *filename_output) 
{
    int i, j;
    FILE *outputFile;

    /* Open Output File */
    outputFile = fopen(filename_output, "w");

    /* Write Output */
    for(i=0;i<nX;i++) {
        for(j=0;j<nY;j++) {
            fprintf(outputFile, "%f ", outputArray[i][j]);
        }
        fprintf(outputFile, "\n");
    }

    /* Close Output File */
    fclose(outputFile);
}
