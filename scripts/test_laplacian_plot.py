#!/usr/bin/env python2
from __future__ import division
import h5py
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot

# A script to plot data saved by the test_laplacian unit test.
# This test computes the Laplacian of a simple function with a FFT method,
# and compares the result to the analytic solution. This script can pinpoint
# errors in the process by showing where the two solutions differ.

# Load data
file = h5py.File("data/test_laplacian.h5")
Nx = file.attrs["grid_sizex"]
Ny = file.attrs["grid_sizey"]
scale = file.attrs["grid_delta"]
Z = file["states"][0,0]  # The computed laplacian is saved in slot 0
Z2 = file["states"][0,1] # The analytic solution is saved in slot 1

X = (2*arange(Nx)-Nx+1)*0.5*scale
Y = (2*arange(Ny)-Ny+1)*0.5*scale
X, Y = meshgrid(X, Y)

# Plot computed Laplacian
fig = pyplot.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet)
# Plot exact Laplacian
fig = pyplot.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, Z2, rstride=1, cstride=1, cmap=cm.bone)
# Plot difference
fig = pyplot.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, abs(Z-Z2), rstride=1, cstride=1, cmap=cm.hot)

pyplot.show()
