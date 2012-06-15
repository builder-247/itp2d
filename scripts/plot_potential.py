#!/usr/bin/env python2
from __future__ import division
from math import *
import sys, os
import h5py
from optparse import OptionParser
from numpy import *

# Plot a stored potential from an itp2d datafile

if __name__=="__main__":
    parser = OptionParser(usage="%prog datafile.h5")
    parser.add_option("-2", "--2d", action="store_true", dest="twodee",
            help="Plot 2D projection of potential")
    parser.add_option("", "--matplotlib3d", action="store_true", dest="matplotlib3d",
            help="Plot 3D image using matplotlib. Very slow.")
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
    if options.twodee or options.matplotlib3d:
        # Plot using matplotlib
        from matplotlib import cm, pyplot, ticker
        from mpl_toolkits.mplot3d import Axes3D
        fig = pyplot.figure(figsize=(12,12))
        ax = fig.add_subplot(1, 1, 1, projection=("rectilinear" if options.twodee else "3d"), aspect="equal")
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
    else:
        # Plot using mayavi
        import mayavi
        from mayavi import mlab
        fig = mlab.figure(bgcolor=(1,1,1))
        light_manager = fig.scene.light_manager
        light_manager.number_of_lights = 4
        lights = light_manager.lights
        for light in lights:
            light.activate = True
            light.intensity = 0.5
        X = transpose(X)
        Y = transpose(Y)
        s = mlab.surf(X, Y, potential, warp_scale='auto')
        mlab.show()
