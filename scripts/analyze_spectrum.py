#!/usr/bin/env python2
# vim: set fileencoding=utf8
from __future__ import division, with_statement
import sys, bisect, math, h5py, numpy
from scipy import *
from scipy.fftpack import fftfreq
from scipy.integrate import trapz, quad, romberg
from numpy.linalg import lstsq
from matplotlib.pyplot import *

# A helper script for computing and plotting nearest-neighbour level spectra,
# spectral rigidity and other statistical measures used in quantum chaos studies.

# The spectral rigidity. Please see, e.g., H.-J. Stöckmann, Quantum chaos: an
# introduction (2006), page 112.
def rigidity(L, e, n):
    if L == 0:
        return 0
    ds = []
    for E in linspace(min(e)+L/2, max(e)-L/2, num=500):
        es = linspace(E-L/2, E+L/2, num=50)
        b, a = polyfit(es, [n(t) for t in es], 1)
        x2s = array([ (n(es[t]) - a - b*es[t])**2 for t in range(len(es)) ])
        d = trapz(x2s, es)/L
        ds.append(d)
    return average(ds)

def main():
    filename = (sys.argv[-1] if sys.argv[-1][-3:] == '.h5' else "data/itp2d.h5")
    file = h5py.File(filename, 'r')
    # read energies
    num_converged = file.attrs["num_converged"]
    energies = array(file["/final_energies"])[:num_converged]
    #energies = range(1000)
    #energies.sort()
    # calculate NND spectrum
    nnd = array([ energies[i] - energies[i-1] for i in range(1, len(energies)) ])
    # normalize average distance to one
    d = mean(nnd)
    nnd /= d
    energies /= d
    # the normalized spectral staircase
    n = lambda t: bisect.bisect(energies, t)

    figure()
    title("Normalized spectral staircase")
    xlabel("$\epsilon$")
    ylabel("$N(\epsilon)$")
    e = linspace(0,max(energies),num=1000)
    plot(e, [ n(t) for t in e ])
    xlim(0, e.max())

    figure()
    title("Spectral rigidity")
    xlabel("$L$")
    ylabel("$\Delta_3(L)$")
    L = linspace(0, 20, num=50)
    D = array([ rigidity(l, e, n) for l in L])
    # Compute rigidity for the picket-fence spectrum
    pftest = array([ rigidity(l, range(1000), lambda t: bisect.bisect(range(1000), t)) for l in L ])
    # Compute rigidity for the Poisson spectrum
    podata = numpy.random.random_sample(size=1000)
    podata.sort()
    ponnd = array([ podata[i] - podata[i-1] for i in range(1, len(podata)) ])
    podata /= mean(ponnd)
    potest = array([ rigidity(l, podata, lambda t: bisect.bisect(podata, t)) for l in L ])
    plot(L, D, label="%d energies from data" % len(energies))
    plot(L, pftest, label="picket-fence (calculated)")
    plot(L, potest, label="Poisson (calculated)")
    plot(L, L/15, ':', label="Poisson")
    plot(L, log(L)/pi**2 - 0.007, ':', label="GOE, $L>>1$")
    plot(L, [1/12]*len(L), ':', label="picket-fence")
    legend(loc="center right")
    xlim(0, max(L))
    ylim(0, 0.4)

    figure()
    title("Nearest neighbour level spacing spectrum")
    xlabel('NNLS')
    ylabel('Frequency')
    n, bins, patches = hist(nnd[nnd <= 6], 300)

    show()

if __name__=="__main__":
    main()
