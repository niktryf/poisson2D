######################################################################
### Python script for visualizing results from poisson2D output
### 
### This script reads grid data from "gridData.txt". 
### The user can simply run and create the 3D plot of the solution by:
### $ python visualize.py
###
### Author: Nikos Tryfonidis
######################################################################
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#################################################################################
### Read Grid Data from "gridData.txt":
#################################################################################
with open('gridData.txt', 'r') as f:
    for line in f:
        strLine = line.split()
        if strLine[0] == 'domain':
            x_0, x_L, y_0, y_L = float(strLine[1]), float(strLine[2]), float(strLine[3]), float(strLine[4])
        if strLine[0] == 'stepsize':
            dx, dy = float(strLine[1]), float(strLine[2])

#################################################################################

# Read whole 2D array into list
with open('output.txt') as f:
    array2d = [[float(digit) for digit in line.split()] for line in f]

# Convert list to numpy array
a = np.asarray(array2d)

# Prepare Surface Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('')

# Set x, y (put grid dimensions and step size here)
X = np.arange(x_0, x_L + dx, dx)
Y = np.arange(y_0, y_L + dy, dy)
# Create meshgrid, set indexing to C-style
X, Y = np.meshgrid(X, Y, indexing='ij')

### Analytic solution for comparison (if known)
#R = np.sin(2*np.pi*(X+Y))
### Print global error (average absolute error) 
#error = np.sqrt((a-R)*(a-R))
#print "Expected order of error: O(%f,%f)" %(dx*dx, dy*dy)
#print "Average Error: %f" %( np.sum(error)/(len(X)*len(Y)) )

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
#############################
ax.set_zlim(np.amin(a), np.amax(a))

plt.show()
