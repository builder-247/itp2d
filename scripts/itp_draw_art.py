#!/usr/bin/env python2
from __future__ import division
import os
import h5py
from math import floor, ceil, sqrt
from numpy import flipud
from optparse import OptionParser
from PIL import Image

# Print pretty pictures from the wave functions computed with itp2d.

colorschemes = {
    "default" : ("RGB", lambda x: (int(256*x), 0, int(256*x*(1-x)))),
    "cold" : ("RGB", lambda x: (int(256*x**2), int(256*x**2), int(256*x**0.75))),
    "finland" : ("RGB", lambda x: (int(256*(1-x)), int(256*(1-x)), 256)),
    "bow" :     ("L",   lambda x: int(256*(1-x)))
}

# Try to import more colorschemes from Matplotlib
try:
    from matplotlib import pylab
    for mapname in pylab.cm.datad:
        colorfunc = lambda x, mapname=mapname: tuple((int(256*t) for t in getattr(pylab.cm, mapname)(x)[:3]))
        colorschemes[mapname] = ("RGB", colorfunc)
except ImportError:
    pass

if __name__ == "__main__":
    # parse command line arguments
    parser = OptionParser(usage="%prog [options] [datafile.h5] [indices]")
    parser.set_defaults(verbose=True, combined=True, colorscheme="default",
            slot=-1, margin=0, trim=0, average_point=0)
    parser.add_option("-v", "--verbose", action="store_true")
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-o", "--output", type="string", metavar="FILENAME", help="Filename for the output image")
    parser.add_option("-a", "--all", action="store_true", help="Draw also the extra states not intended to converge")
    parser.add_option("-s", "--square", action="store_true", help="Discard states so that the resulting image is square")
    parser.add_option("-m", "--margin", type="int", metavar="PIXELS", help="Add some empty space between states in the combined image")
    parser.add_option("-t", "--trim", type="int", metavar="PIXELS", help="Remove some space from the edges of each state image")
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
    assert (N == states.shape[1]), "Datafile corrupted, recorded num_states does not match state array dimensions."
    if ((options.slot >= 0) and options.slot >= states.shape[0]) \
            or ((options.slot < 0) and abs(options.slot) > states.shape[0]):
        parser.error("Invalid slot index %d. Datafile has %d slots." % (options.slot, states.shape[0]))
    trim = options.trim
    Mx = file.attrs["grid_sizex"] - 2*trim
    My = file.attrs["grid_sizey"] - 2*trim
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
    if options.combined:
        full_im = Image.new(mode, (columns*Mx+(columns-1)*margin, rows*My+(rows-1)*margin), colorfunc(0))
    state_im = Image.new(mode, (Mx, My), colorfunc(0))
    counter = 0
    # Loop through all states to be plotted
    for index in indices:
        # The data to plot is the square of the absolute value of the
        # wave function, i.e., the probability density
        Z = abs(states[options.slot, index][trim:-trim,trim:-trim])**2
        # Normalize
        if options.average_point == 0:
            Z /= Z.max()
        else:
            avg = Z.mean()
            Z *= options.average_point/avg
        # Flatten the data array, map each value to a color, and write to the image
        # It seems that PIL writes data in reverse order (bottom up), so we need to
        # flip the data array with flipud
        state_im.putdata([ colorfunc(x) for x in flipud(Z).flat ])
        if options.combined:
            paste_corner = ((counter % columns)*(Mx+margin), (counter // columns)*(My+margin))
            full_im.paste(state_im, paste_corner)
        if options.separate:
            out = out_multifilename % index
            state_im.save(out)
            if options.verbose:
                print "State %d drawn. Saved in %s" % (index, out)
        elif options.verbose:
            print "State %d drawn." % index
        counter += 1
    if options.combined:
        full_im.save(out_filename)
        if options.verbose:
            print "Combined picture saved in %s" % out_filename
