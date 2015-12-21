######################################################################
### Python script for visualizing results from poisson2D output
### 
### User must enter 2 command line arguments, the step sizes 
### for x and y respectively. These can be found in the gridData.txt
### text file:
### 
### $ python visualize.py <dx> <dy>
###
### These are used to calculate the axis scales.
###
### Author: Nikos Tryfonidis
######################################################################
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


#################################################################################
### Command line arguments:
#################################################################################

# Check number of command line arguments
if len(sys.argv) != 3:
    print "Usage: python visualize.py <dx> <dy>"
    print "Values for dx and dy can be found in file 'gridData.txt' "
    print "Please run again following the command line input format above."
    print "Exiting..."
    sys.exit(1)

# Get command line arguments
dx = float(sys.argv[1])
dy = float(sys.argv[2])
#################################################################################


# Read whole 2D array into list
with open('output.txt') as file:
    array2d = [[float(digit) for digit in line.split()] for line in file]

# Convert list to numpy array
a = np.asarray(array2d)

# Prepare Surface Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Numerical Solution')

# Set x, y (put grid dimensions and step size here)
X = np.arange(0, 1 + dx, dx)
Y = np.arange(0, 1 + dy, dy)
X, Y = np.meshgrid(X, Y)

# Analytic solution for comparison (if known)
#R = np.sin(2*np.pi*(X+Y))

# Create surface plot
# rstride, cstride: the input stride for rows and columns respectively,
# useful for large grids.

surface = ax.plot_surface(X, Y, a, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$u$')

#############################
# Set z axis range
ax.set_zlim(np.amin(a), np.amax(a))

plt.show()
