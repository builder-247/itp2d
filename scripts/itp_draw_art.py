#!/usr/bin/env python2
from __future__ import division
import os
from math import floor, ceil, sqrt
from collections import namedtuple
import numpy as np
import h5py
from optparse import OptionParser
from PIL import Image, ImageFont, ImageDraw

# Print pretty pictures from the wave functions computed with itp2d.

try:
    from progressbar import ProgressBar, ETA, Percentage, Bar
    has_progressbar = True
except ImportError:
    # Provide mock objects in case import fails
    has_progressbar = False
    print "Import of module 'progressbar' failed. Will not print a progress bar."
    class ProgressBar(object):
        def __init__(self, *args, **kwargs):
            pass
        def __call__(self, x):
            return x
    class ETA(object):
        pass
    class Percentage(object):
        pass
    class Bar(object):
        pass

# Each color scheme is a tuple with two attributes: a mode (either "RGB" or
# "L") and a function. The function maps a NxN numpy array with float values to
# a NxNx3 numpy array of uint8 RGB values (if mode is "RGB") or to a NxN numpy
# array of uint8 values (if mode is "L"). The values of the original array are
# assumed to be in the range [0, 1]. Everything beyond this range is clipped.
Colorscheme = namedtuple("Colorscheme", ["mode", "func"])

CC = conv_const = np.nextafter(256,0) # Conversion constant from float to uint8

def C(x):
    return np.clip(x, 0, 1)

# Sorry for the ugly function definitions... This asks for some Haskell-style
# function composition
colorschemes = {
    "default" : Colorscheme("RGB",
        lambda x: np.uint8(CC*np.dstack((C(x), np.zeros_like(C(x)), C(x)*(1-C(x)))))),
    "cold" : Colorscheme("RGB",
        lambda x: np.uint8(CC*np.dstack((C(x)**2, C(x)**2, C(x)**0.75)))),
    "finland" : Colorscheme("RGB",
        lambda x: np.uint8(CC*np.dstack((1-C(x), 1-C(x), np.ones_like(C(x)))))),
    "bow" : Colorscheme("L",
        lambda x: np.uint8(CC*(1-C(x))))
}

# Try to import more colorschemes from Matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pylab
    for mapname in pylab.cm.datad:
        colorfunc = lambda x, cmap=getattr(pylab.cm, mapname): cmap(x, bytes=True)[...,:3]
        colorschemes[mapname] = Colorscheme("RGB", colorfunc)
except ImportError:
    pass

if __name__ == "__main__":
    # parse command line arguments
    parser = OptionParser(usage="%prog [options] [datafile.h5] [indices]")
    parser.set_defaults(verbose=True, labels=False, label_font_size=10,
            combined=True, colorscheme="default", slot=-1, margin=0, trim=0,
            average_point=0, rescale=1)
    parser.add_option("-v", "--verbose", action="store_true")
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-o", "--output", type="string", metavar="FILENAME", help="Filename for the output image")
    parser.add_option("-a", "--all", action="store_true", help="Draw also the extra states not intended to converge")
    parser.add_option("-l", "--labels", action="store_true", help="Draw labels with each state's index and energy to the combined image")
    parser.add_option("", "--label-font-size", type="int", help="Font size for the labels")
    parser.add_option("-s", "--square", action="store_true", help="Discard states so that the resulting image is square")
    parser.add_option("-r", "--rescale", type="float", metavar="FACTOR", help="Rescale resulting images with this factor")
    parser.add_option("-m", "--margin", type="int", metavar="PIXELS", help="Add some empty space between states in the combined image (after rescaling)")
    parser.add_option("-t", "--trim", type="int", metavar="PIXELS", help="Remove some space from the edges of each state image (before rescaling)")
    parser.add_option(      "--separate", action="store_true", help="Save separate images of each state")
    parser.add_option(      "--no-combined", action="store_false", dest="combined", help="Do not save a combined image of all the states")
    parser.add_option("-c", "--colorscheme", type="choice", choices=colorschemes.keys(), help="The colors. Valid choices: %s" % colorschemes.keys())
    parser.add_option("-S", "--slot", type="int", help="If the datafile contains several sets of states this lets you select the one you want. The default is to load the last one.")
    parser.add_option(      "--scale-to-average-point", type="float", metavar="VAL", dest="average_point", default=0,
            help="Instead of scaling density data so that 1.0 corresponds to maximum density, scale so that VAL corresponds to average density.")
    (options, args) = parser.parse_args()
    if (len(args) > 0):
        filename = args[0]
    else:
        filename = "data/itp2d.h5"
    indices = []
    for arg in args[1:]:
        try:
            index = int(arg)
        except ValueError:
            parser.error("Could not interpret '%s' as an index (should be a number)" % arg)
        if (index < 0):
            parser.error("Could not interpret '%s' as an index (should be a non-negative)" % arg)
        indices.append(index)
    if options.output:
        out_filename = options.output
    else:
        out_filename = os.path.splitext(filename)[0]+".png"
    out_multifilename = os.path.splitext(out_filename)[0]+"%d.png"
    margin = options.margin
    # Option parsing done. Read data.
    file = h5py.File(filename, 'r')
    N = file.attrs["num_states"]
    states = file["/states"]
    energies = file["/final_energies"]
    assert (N == states.shape[1]), "Datafile corrupted, recorded num_states does not match state array dimensions."
    if ((options.slot >= 0) and options.slot >= states.shape[0]) \
            or ((options.slot < 0) and abs(options.slot) > states.shape[0]):
        parser.error("Invalid slot index %d. Datafile has %d slots." % (options.slot, states.shape[0]))
    trim = options.trim
    Mx = file.attrs["grid_sizex"] - 2*trim
    My = file.attrs["grid_sizey"] - 2*trim
    scaled_Mx = int(options.rescale * Mx)
    scaled_My = int(options.rescale * My)
    # The states that are drawn are either
    # * The ones specified by the indices argument
    # * All states found in the file (if --all)
    # * Only the ones required to converge
    # In addition, if --square is set, cut the number of states so that
    # we get a nice square image
    if len(indices) > 0:
        num_to_draw = len(indices)
    elif options.all:
        num_to_draw = N
        indices = range(N)
    else:
        num_to_draw = file.attrs["num_wanted_to_converge"] + file.attrs["ignore_lowest"]
        indices = range(num_to_draw)
    if options.square:
        num_to_draw = int(sqrt(num_to_draw))**2
        indices = indices[:num_to_draw]
    columns = int(floor(sqrt(num_to_draw)))
    rows = int(ceil(num_to_draw/columns))
    mode, colorfunc = colorschemes[options.colorscheme]
    # Load PIL's default font for creating labels
    font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSansMono.ttf",
            options.label_font_size)
    # Get color of text from the colormap
    if mode == "RGB":
        foreground_color = tuple(colorfunc(np.ones((1,1)))[0,0])
    else:
        foreground_color = colorfunc(np.ones((1,1)))[0,0]
    # Initialize combined image
    if options.combined:
        # Get color of zero from colormap
        if mode == "RGB":
            background_color = tuple(colorfunc(np.zeros((1,1)))[0,0])
        else:
            background_color = colorfunc(np.zeros((1,1)))[0,0]
        full_x = columns*scaled_Mx+(columns-1)*margin
        full_y = rows*scaled_My+(rows-1)*margin
        full_im = Image.new(mode, (full_x, full_y), background_color)
        full_draw = ImageDraw.Draw(full_im)
    counter = 0
    # Loop through all states to be plotted
    progressbar = ProgressBar(widgets=["Drawing images: ", Percentage(), Bar(), ETA()])
    for index in progressbar(indices):
        # The data to plot is the square of the absolute value of the
        # wave function, i.e., the probability density
        if options.trim == 0:
            Z = abs(states[options.slot, index])**2
        else:
            Z = abs(states[options.slot, index][trim:-trim,trim:-trim])**2
        # Flip, since in the original data array y-axis points "downwards"
        Z = np.flipud(Z)
        # Normalize
        if options.average_point == 0:
            Z /= Z.max()
        else:
            avg = Z.mean()
            Z *= options.average_point/avg
        # Create image by mapping data through colorfunc
        state_im = Image.fromarray(colorfunc(Z), mode=mode)
        if options.rescale != 1:
            state_im.thumbnail((scaled_Mx, scaled_My), Image.ANTIALIAS)
        if options.combined:
            paste_x = (counter % columns)*(scaled_Mx+margin)
            paste_y = (counter // columns)*(scaled_My+margin)
            paste_corner = (paste_x, paste_y)
            full_im.paste(state_im, paste_corner)
            if options.labels:
                E = energies[index]
                label = "n = %d, E = %.3f" % (index, E)
                full_draw.text(paste_corner, label, fill=foreground_color, font=font)
        if options.separate:
            out = out_multifilename % index
            state_im.save(out)
            if options.verbose and not has_progressbar:
                print "State %d drawn. Saved in %s" % (index, out)
        elif options.verbose and not has_progressbar:
            print "State %d drawn." % index
        counter += 1
    if options.combined:
        full_im.save(out_filename)
        if options.verbose:
            print "Combined picture saved in %s" % out_filename
