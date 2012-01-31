#!/usr/bin/env python
from __future__ import division, with_statement
import sys, os
import h5py
from numpy import *
from optparse import OptionParser
from PIL import Image, ImageFilter

# Print pretty pictures from the wave functions computed with itp2d.

colorschemes = {
    "default" : ("RGB", lambda x: (int(256*x), 0, int(256*x*(1-x)))),
    "finland" : ("RGB", lambda x: (int(256*x**2), int(256*x**2), int(256*x**0.75))),
    "bow" :     ("L",   lambda x: int(256*(1-x)))
}

if __name__ == "__main__":
    # parse command line arguments
    parser = OptionParser(usage="%prog [options] [datafile.h5] [indices]")
    parser.set_defaults(verbose=True, combined=True, colorscheme="default", slot=-1)
    parser.add_option("-v", "--verbose", action="store_true")
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-o", "--output", type="string", metavar="FILENAME", help="Filename for the output image")
    parser.add_option("-a", "--all", action="store_true", help="Draw also the extra states not intended to converge")
    parser.add_option("-s", "--square", action="store_true", help="Discard states so that the resulting image is square")
    parser.add_option(      "--separate", action="store_true", help="Save separate images of each state")
    parser.add_option(      "--no-combined", action="store_false", dest="combined", help="Do not save a combined image of all the states")
    parser.add_option("-c", "--colorscheme", type="choice", choices=colorschemes.keys(), help="The colors. Valid choices: %s" % colorschemes.keys())
    parser.add_option("-S", "--slot", type="int", help="If the datafile contains several sets of states this lets you select the one you want. The default is to load the last one.")
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
    # Option parsing done. Read data.
    file = h5py.File(filename, 'r')
    N = file.attrs["num_states"]
    states = file["/states"]
    assert (N == states.shape[1]), "Datafile corrupted, recorded num_states does not match state array dimensions."
    if ((options.slot >= 0) and options.slot >= states.shape[0]) \
            or ((options.slot < 0) and abs(options.slot) > states.shape[0]):
        parser.error("Invalid slot index %d. Datafile has %d slots." % (options.slot, states.shape[0]))
    Mx = file.attrs["grid_sizex"]
    My = file.attrs["grid_sizey"]
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
        full_im = Image.new(mode, (columns*Mx, rows*My))
    state_im = Image.new(mode, (Mx, My))
    counter = 0
    # Loop through all states to be plotted
    for index in indices:
        # The data to plot is the square of the absolute value of the
        # wave function, i.e., the probability density
        Z = abs(states[options.slot, index])**2
        # Normalize
        Z /= Z.max()
        # Flatten the data array, map each value to a color, and write to the image
        state_im.putdata([ colorfunc(x) for x in Z.flat ])
        if options.combined:
            paste_corner = ((counter % columns)*Mx, (counter // columns)*My)
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
