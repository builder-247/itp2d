#!/usr/bin/env python
from numpy import *
import matplotlib.pyplot as pyplot
import matplotlib.mlab as mlab
import h5py

# Plot data produced by the unit tests testing the random
# number generator. The unit tests produce normally and
# uniformly distributed random numbers and tests that the
# mean and standard deviation calculated from these numbers
# are correct. This script can be used to visually confirm
# what the distributions look like.

pyplot.figure()
datafile = h5py.File("data/test_rng_gaussian.h5")
x = datafile["rand"]
n, bins, patches = pyplot.hist(x, 100, normed=1)
y = mlab.normpdf(bins, 0.0, 1.0)
pyplot.plot(bins, y, 'r', linewidth=1)

pyplot.figure()
datafile = h5py.File("data/test_rng_uniform.h5")
x = datafile["rand"]
n, bins, patches = pyplot.hist(x, 100, normed=1)
y = [1.0]*len(bins)
pyplot.plot(bins, y, 'r', linewidth=1)

pyplot.show()
