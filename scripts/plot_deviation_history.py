#!/usr/bin/env python2
from __future__ import division
import sys, os, h5py
from optparse import OptionParser
import matplotlib
from matplotlib import pyplot
from numpy import *

# Plots how the energy standard deviation or other measures of accuracy
# evolve during iterations of ITP.

NUM = 12
nan = float("nan")

# A simple enum for what to plot
class Mode:
    relative = object()
    absolute = object()
    final = object()
    change = object()

def main():
    # Parse command line arguments
    parser = OptionParser(usage="%prog [options] datafile.h5")
    parser.set_defaults(mode=Mode.relative, energy_relative=False)
    parser.add_option("-r", "--energy_relative", action="store_const", const=Mode.relative, dest="mode",
            help="treat standard deviations relative to energy")
    parser.add_option("-a", "--absolute", action="store_const", const=Mode.absolute, dest="mode",
            help="plot absolute standard deviations instead of relative ones")
    parser.add_option("-f", "--final", action="store_const", const=Mode.final, dest="mode",
            help="plot difference to final value at each step")
    parser.add_option("-c", "--change", action="store_const", const=Mode.change, dest="mode",
            help="plot relative changes between successive steps")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("required argument datafile missing")
    filename = args[0]
    file = h5py.File(filename)
    energyhist = array(file["/energy_history"])
    hist = array(file["/deviation_history"])
    try:
        # import time_step history if available
        epshist = array(file["/time_step_history"])
    except:
        pass
    file.close()
    steps = range(1,hist.shape[0]+1)
    # Choose some some representative states to follow; do not plot all of them
    states_to_plot = [ int(x) for x in linspace(0, hist.shape[1]-1, num=NUM) ]
    for n in states_to_plot:
        if options.mode == Mode.relative:
            pyplot.ylabel("relative $\sigma(E)$ (a.u.)")
            dE = [ hist[i,n]/energyhist[i,n] for i in range(hist.shape[0]) ]
            pyplot.semilogy(steps[1:], dE[1:], label=str(n))
        elif options.mode == Mode.absolute:
            pyplot.ylabel("absolute $\sigma(E)$ (a.u.)")
            dE = [ hist[i,n] for i in range(hist.shape[0]) ]
            pyplot.semilogy(steps[1:], dE[1:], label=str(n))
        elif options.mode == Mode.final:
            pyplot.ylabel("$|\sigma(E)-\sigma(E)_f|$ (a.u.)")
            dE = [ abs(hist[i,n] - hist[-1,n]) for i in range(hist.shape[0]) ]
            pyplot.semilogy(steps, dE, label=str(n))
        elif options.mode == Mode.change:
            pyplot.ylabel("relative $|\Delta\sigma(E)|$ (a.u.)")
            dE = [nan]+[ abs((hist[i,n]-hist[i-1,n])/hist[i-1,n]) for i in range(1,hist.shape[0]) ]
            pyplot.plot(steps, dE, label=str(n))
        else:
            raise RuntimeError("options.mode had unexpected value")
    # Plot also values of the time step if that data is available
    if 'epshist' in dir():
        for (step, eps) in epshist:
            if (step > hist.shape[0]):
                continue
            pyplot.axvline(step, linewidth=1, linestyle=':', color='k')
            pyplot.text(step+0.2, pyplot.gca().get_ylim()[1], "%.1e" % eps, verticalalignment='bottom', size='small')
    pyplot.xlim(1, max(steps))
    pyplot.xlabel("Step")
    pyplot.legend(loc="lower left", ncol=4, prop={"size": "small"})
    pyplot.show()

if __name__=="__main__":
    main()
