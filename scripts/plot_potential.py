#!/usr/bin/env python2
from __future__ import division
from math import *
import sys, os
import h5py
from optparse import OptionParser
from numpy import *

try:
    from matplotlib import cm, pyplot, ticker
    from mpl_toolkits.mplot3d import Axes3D
except ImportError as e:
    print "Couldn't import matplotlib:", str(e)
    has_matplotlib = False
else:
    has_matplotlib = True

try:
    import mayavi
    from mayavi import mlab
except ImportError as e:
    print "Couldn't import mayavi:", str(e)
    has_mayavi = False
else:
    has_mayavi = True

# Plot a stored potential from an itp2d datafile
def main():
    parser = OptionParser(usage="%prog datafile.h5")
    parser.add_option("", "--clip", type="float", help="Clip potential values larger than this")
    parser.add_option("-o", "--output", type="string", metavar="FILE",
            help="Save image to file instead of showing")
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
    mask = (None if options.clip == None else (potential > options.clip))

    X = (2*arange(Mx)-Mx+1)*0.5*scale
    Y = (2*arange(My)-My+1)*0.5*scale
    X1d = X[:]
    Y1d = Y[:]
    extent = array([min(X), max(X), min(Y), max(Y)])
    X, Y = meshgrid(X, Y)
    # Actual plotting
    if options.twodee or options.matplotlib3d:
        if options.clip is not None:
            clip(potential, float("-inf"), options.clip, out=potential)
        matplotlib_plot(X, Y, potential, options, extent)
    else:
        mayavi_3dplot(X1d, Y1d, transpose(potential), options, transpose(mask))

def matplotlib_plot(X, Y, Z, options, extent):
        fig = pyplot.figure(figsize=(12,12))
        ax = fig.add_subplot(1, 1, 1,
                projection=("rectilinear" if options.twodee else "3d"), aspect="equal")
        if options.twodee:
            ax.imshow(Z, extent=extent, rasterized=True, origin='lower')
        else:
            ax.set_xlabel('$x$')
            ax.set_ylabel('$y$')
            ax.plot_surface(X, Y, Z,
                    rstride=1, cstride=1, cmap=cm.jet, linewidth=0)
            ax.apply_aspect()
        if options.output:
            pyplot.savefig(options.output)
        else:
            pyplot.show()

def mayavi_3dplot(X, Y, Z, options, mask):
    if not has_mayavi:
        raise RuntimeError("Can't plot with mayavi because it is not available")
    fig = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))
    light_manager = fig.scene.light_manager
    light_manager.number_of_lights = 4
    lights = light_manager.lights
    for light in lights:
        light.activate = True
        light.intensity = 0.5
    s = mlab.surf(X, Y, Z, warp_scale='auto', mask=mask)
    #mlab.axes(s, nb_labels=10)
    if options.output:
        mlab.savefig(options.output)
    else:
        mlab.show()
    mlab.close(all=True)

if __name__=="__main__":
    main()
