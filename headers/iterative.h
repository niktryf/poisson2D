double **jacobiIteration2D(double **u, double **f, 
                           double dx, double dy, int nX, int nY);

double **weightedJacobiIteration2D(double **u, double **f, 
                           double dx, double dy, int nX, int nY);

double **gaussSeidelIteration2D(double **u, double **f, 
                                double dx, double dy, int nX, int nY);

double **redBlackGaussSeidelIteration2D(double **u, double **f, 
                                double dx, double dy, int nX, int nY);

double **SORIteration2D(double **u, double **f, 
                        double dx, double dy, int nX, int nY);

double **redBlackSORIteration2D(double **u, double **f, 
                           double dx, double dy, int nX, int nY);
