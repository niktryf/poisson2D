***************************************************************************            
            Poisson2D: A 2D Poisson solver in C
            Author: Nikos Tryfonidis
            
            The MIT License (MIT)
            Copyright (c) 2015 Nikos Tryfonidis
***************************************************************************

This project is a numerical solver for the Poisson PDE in 
two dimensions. The current README contains brief information
on compiling and running.

The code has been parallelized with OpenMP and is maintained with the 
OpenMP additions.

For further information, please see the documentation in the 
"documentation" directory.

---------------------------------------------------------------------------
Compilation

The program should compile without problems on any linux distribution 
with a C compiler. To compile, simply use the provided makefile, by typing 
"make" in the main directory, where the makefile is located.

To remove the executable, type "make clean".

----------------------------------------------------------------------------
Execution

The executable, named poisson2D, needs two arguments: 
1. Number of grid points in X
2. Number of grid points in Y

For example, to run with 101 grid points in both dimensions, run:
./poisson2D 101 101

It is recommended to specify the desired number of OpenMP threads before running,
for example:

export OMP_NUM_THREADS=2

for two threads.

In case your computer uses hyperthreading, it is recommended to tie threads to 
different physical cores, to keep the system from using hyperthreaded cores. 
This can be done by using, for example:

export GOMP_CPU_AFFINITY="0 1 2 3"

The above assigns processes 0 to 3 to different threads. Assuming that these are 
different physical cores, this should avoid using hyperthreading, which performs 
poorly for applications as this.

---------------------------------------------------------------------------------
Output

Output is written in the "output" directory. The solution is in "output.txt".
"gridData.txt" contains grid data read by the plotting python script 
"visualize.py".

You can produce a 3D surface plot of the solution by running the python 
script. 

Simply type:

python visualize.py

in the output directory.

Please note that the python script needs numpy and matplotlib, two popular 
numerical packages for Python.

For further details, please see the documentation.





