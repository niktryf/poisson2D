double **solvePoisson2D(double ** (*iterate)(double **, double **, double, double, int, int), 
                        double **array2D, double **rhs, double dx, double dy, int nX, int nY);
