#!/usr/bin/env python2
# vim: set fileencoding=utf8
from __future__ import division
import h5py
import warnings
from argparse import ArgumentParser
import numpy as np

try:
    from matplotlib import pyplot
except ImportError:
    warnings.warn("Could not import matplotlib, 2D plot not available")

try:
    from mayavi import mlab
except ImportError:
    warnings.warn("Could not import mayavi, 3D plot not available")

# Plot a stored potential from an itp2d datafile
def main():
    parser = ArgumentParser()
    parser.add_argument("-c", "--clip", type=float, help="Clip potential values larger than this")
    parser.add_argument("-o", "--output", type=str, metavar="FILE",
            help="Save image to file instead of showing")
    parser.add_argument("-2", "--2d", action="store_true", dest="twodee",
            help="Plot 2D projection of potential with matplotlib")
    parser.add_argument("--clip-hack-number-of-contours", type=int, default=100,
            help="How many contours to use as the surface to give nice clipping in mayavi")
    parser.add_argument("datafile", type=str, help="Datafile saved by itp2d")
    args = parser.parse_args()

    # Load data
    with h5py.File(args.datafile, 'r') as f:
        sx = f.attrs["grid_sizex"]
        sy = f.attrs["grid_sizey"]
        grid_delta = f.attrs["grid_delta"]
        potential = f["potential_values"].value

    X = (2*np.arange(sx)-sx+1)*0.5*grid_delta
    Y = (2*np.arange(sy)-sx+1)*0.5*grid_delta
    extent = np.array([min(X), max(X), min(Y), max(Y)])
    # Actual plotting
    if args.twodee:
        if args.clip is not None:
            np.clip(potential, float("-inf"), args.clip, out=potential)
        matplotlib_2dplot(potential, args, extent)
    else:
        mayavi_3dplot(X, Y, np.transpose(potential), args)

def matplotlib_2dplot(Z, args, extent):
    if not "pyplot" in globals():
        raise RuntimeError("Can't plot with matplotlib because it's not available")
    fig = pyplot.figure(figsize=(12,12))
    ax = fig.add_subplot(1, 1, 1, projection="rectilinear")
    ax.imshow(Z, extent=extent, rasterized=True, origin='lower')
    if args.output:
        pyplot.savefig(args.output)
    else:
        pyplot.show()

def mayavi_3dplot(X, Y, Z, args):
    if not "mlab" in globals():
        raise RuntimeError("Can't plot with mayavi because it is not available")
    fig = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))
    light_manager = fig.scene.light_manager
    light_manager.number_of_lights = 4
    lights = light_manager.lights
    for light in lights:
        light.activate = True
        light.intensity = 0.5
    s = mlab.surf(X, Y, Z, warp_scale='auto')
    if args.clip is not None:
        s.module_manager.scalar_lut_manager.reverse_lut = True
        s.module_manager.scalar_lut_manager.use_default_range = False
        s.module_manager.scalar_lut_manager.data_range = np.array([0,  args.clip])
        s.contour.number_of_contours = args.clip_hack_number_of_contours
        s.enable_contours = True
        s.contour.filled_contours = True
        s.contour.maximum_contour = args.clip
    if args.output:
        mlab.savefig(args.output)
    else:
        mlab.show()
    mlab.close(all=True)

if __name__=="__main__":
    main()
