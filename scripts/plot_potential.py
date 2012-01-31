#!/usr/bin/env python
from __future__ import division
from math import *
import sys, os
import h5py
from optparse import OptionParser
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, pyplot, ticker
from numpy import *

# Plot a stored potential from an itp2d datafile

if __name__=="__main__":
    parser = OptionParser(usage="%prog datafile.h5")
    parser.add_option("-2", "--2d", action="store_true", dest="twodee")
    (options, args) = parser.parse_args()

    try:
        filename = args[0]
    except:
        parser.print_usage()
        sys.exit()

    # Option parsing done
    file = h5py.File(filename, 'r')

    Mx = file.attrs["grid_sizex"]
    My = file.attrs["grid_sizey"]
    scale = file.attrs["grid_delta"]
    potential = file["potential_values"].value

    X = (2*arange(Mx)-Mx+1)*0.5*scale
    Y = (2*arange(My)-My+1)*0.5*scale
    X1d = X[:]
    Y1d = Y[:]
    extent = array([min(X), max(X), min(Y), max(Y)])
    X, Y = meshgrid(X, Y)
    # Actual plotting
    fig = pyplot.figure(figsize=(12,12))
    if options.twodee:
        projection = "rectilinear"
    else:
        projection = "3d"
    ax = fig.add_subplot(1, 1, 1, projection=projection, aspect="equal")
    if options.twodee:
        ax.xaxis.set_major_locator(ticker.NullLocator())
        ax.yaxis.set_major_locator(ticker.NullLocator())
        ax.imshow(potential, extent=extent, rasterized=True)
    else:
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.plot_surface(X, Y, potential, rstride=1, cstride=1, cmap=cm.jet, linewidth=0)
        ax.apply_aspect()
pyplot.show()
