itp2d – the fast Schrödinger equation solver for 2D systems
=============================================================

Introduction
------------

itp2d is a fast Schrödinger equation solver for single-particle systems in
two dimensions. It uses the imaginary time propagation (ITP) algorithm,
combined with any-order operator factorization and exact, gauge-invariant
factorization for magnetic fields. The program has been designed to be general,
in the sense that it does not rely on special symmetries of the system or
special basis expansions. It is also tuned for efficiency, especially for
computing the eigenstates and eigenenergies of a system up to very highly
excited states. The structure of the program is made simple and maintainable by
an object-oriented design, implemented in C++. Computations are made fast by
offloading as much work as possible to heavily optimized basic linear algebra
subprograms via standard [CBLAS][] and [LAPACK][] interfaces, and to [FFTW][],
the fastest discrete Fourier transform library in the West.

A user intending to use itp2d for science should first read the [associated
article][article] published in [Computer Physics Communications][CPC]. The
article describes the underlying algorithms behind the program in more detail
and from a scientific point of view, whereas the purpose of this document is to
give a quick overview about the structure and use of the program.

If you use itp2d for published scientific work, please cite the
aforementioned [article][]. You are not bound by any license or user agreement
to do so – this is only a kind wish from a scientist to another. We all have
indices to increase and CVs to polish.

[CBLAS]: http://netlib.org/blas
[LAPACK]: http://www.netlib.org/lapack
[FFTW]: http://www.fftw.org
[CPC]: http://www.journals.elsevier.com/computer-physics-communications
[article]: http://TODO_FILL_IN_THIS_LINK_ONCE_THE_ARTICLE_IS_PUBLISHED

Acquiring itp2d
-----------------

The easiest way to get up-to-date versions of itp2d is to use [GitHub][],
which is a site built for distributing software using the amazing version
control system [Git][]. By using itp2d's [GitHub site][webpage] you can see
recent changes made to the program, report and track bugs found in the program,
access user-generated documentation and even create your own versions of
itp2d.

[github]: https://github.com
[git]: http://git-scm.com

Program license
---------------

itp2d is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

itp2d is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
itp2d.  If not, see <http://www.gnu.org/licenses/>.

[author]: mailto:perttu.luukko@iki.fi
[webpage]: https://github.com/nuteater/itp2d

Program structure
-----------------

### Building the program

#### Requirements

itp2d is only tested on Linux computers and on x86 and x86-64 architectures, but it should
work in any POSIX compliant, little endian system with the software specified below. Porting itp2d
to non-POSIX systems should be easy (should the need ever arise), since itp2d includes very little
platform-dependent code.

In order to build itp2d you must have:

- A fairly recent C++ compiler. The program is written in standards-compliant
  C++, but itp2d is only tested with the [GNU Compiler Collection][gcc]. The
  program also uses some [TR1][] additions to C++, so your C++ standard library
  (which is usually bundled with the compiler) must have a sufficiently
  complete implementation of [TR1][]. If you wish to make use of several
  processors the compiler must also support [OpenMP][]. For GCC, versions 4.3
  and above seem to fulfill these requirements.
- The [Templatized C++ Command Line Parser Library][TCLAP].
- [FFTW][] 3.x, which is used for computing discrete Fourier transforms.
- The [HDF5][] library, version 1.8.x, which itp2d uses for saving data on
  disk. The library must be compiled with C++ support.
- A linear algebra library with a [CBLAS][] and [LAPACK][] interfaces. Possible alternatives are,
  e.g., [ATLAS][], Intel's [MKL][] or AMD's [ACML][]
- [Python][], which is used in generating dependencies automatically during the build.
- To build the unit test functions distributed with itp2d you'll also need
  the [Google C++ Testing Framework (gtest)][gtest]. There is a helper script
  distributed with itp2d which can take care of the installation of gtest.

itp2d also comes with several helper scripts for data analysis and plotting. These scripts can be
found in the `scripts` subfolder. The scripts are written in [Python][], and they require various
Python libraries to function. However, these scripts are not necessary for the normal operation of
itp2d.

[OpenMP]: http://openmp.org
[TCLAP]: http://tclap.sourceforge.net
[gcc]: http://gcc.gnu.org
[hdf5]: http://www.hdfgroup.org/HDF5
[ATLAS]: http://math-atlas.sourceforge.net
[MKL]: http://software.intel.com/en-us/articles/intel-mkl
[ACML]: http://developer.amd.com/libraries/acml/pages/default.aspx
[Python]: http://python.org
[gtest]: http://code.google.com/p/googletest
[TR1]: http://en.wikipedia.org/wiki/C%2B%2B_Technical_Report_1

#### Compiling itp2d

The program comes with a GNU Makefile, so you can compile itp2d simply by
running the command `make` in the directory where the Makefile resides. The
Makefile assumes you are using the `g++` compiler from the GNU Compiler
Collection. If you wish to compile itp2d with some other compiler, you need
to adjust some command line flags used in the Makefile. Please note that the
compiler used should not make a dramatic difference in the performance of itp2d,
since most heavy lifting is delegated to external routines.

You can also build and run the unit test suite by running `make check` and
convert this README file to HTML by running `make doc`. For compiling the unit
tests you need to have the source code of gtest. [The developers of gtest do
not recommend installing gtest as a library][gtestnote] – the gtest source code
should be compiled for all projects using gtest separately. The Makefile
supplied with itp2d assumes you have the gtest source code unpacked to
directory `src/gtest`. Instead of manually unpacking, you can use the provided
script `scripts/fetch_gtest.sh`.

[gtestnote]: http://groups.google.com/group/googletestframework/browse_thread/thread/668eff1cebf5309d

If you have trouble compiling itp2d, for example because one of the required
libraries is not found even though you have it installed, please look into the
Makefile and modify the compiler flags if needed. The Makefile respects your
choice of command line arguments to the C++ compiler, as specified by the
`CXXFLAGS` environment variable.

##### Using MKL or ACML for linear algebra

By default itp2d uses [CBLAS][] and [LAPACK][] for linear algebra. If you wish to use Intel's
[MKL][], edit the Makefile, or create a file with the name `local.mk` in the same directory as the
Makefile with the following content:

	lalib := MKL

Similarly, a file with

	lalib := ACML

instructs the Makefile to compile against AMD's [ACML][] library. You need to recompile itp2d to
this change to take effect.

### Command line parameters

Please run `itp2d --help` to access the embedded documentation about the possible command line
parameters you can pass to itp2d.

### Helper scripts

There are several helper scripts distributed with itp2d, and these scripts you can find in the
`scripts` subdirectory. Most of these scripts are for a very specific purpose, and they are not intended
to be a complete data analysis package. However, they can be a starting point for creating
your own data analysis tools based on itp2d. If you come up with useful tools, please add them as
a patch via the [program webpage][webpage] or send them to the program [author][], and
they will be added to the main itp2d distribution.

The scripts are not required to access the data created by itp2d, since the datafiles created by
itp2d are standard [HDF5][] files with a simple layout. You should be able to access them
easily with any program cabable of reading [HDF5][] files, such as Python or MATLAB.

Examples & use cases
--------------------

### Getting help and reporting bugs

If you have any questions about itp2d, please do not hesitate to contact the [author][]. For
reporting bugs the [GitHub][] issue interface found on the itp2d main [webpage][] is preferable.

### Running unit tests

To build and run the unit tests distributed with itp2d run `make check`. This
creates and runs a binary `run_tests`, This program runs several tests to make
sure everything in the itp2d codebase is functioning properly. You should
always run the tests after making any changes in the itp2d code, and before
running any important simulations!

### Solving the harmonic oscillator

The harmonic oscillator is a nice example of a system with an easy and well known energy spectrum
and eigenstates. When the potential is (in atomic units)

	V(r) = ½r²,

where `r` is the distance from the origin, the energies of the system will be integers. The ground
state energy is 1, the second state is doubly degenerate with energy 2, the third is triply
degenerate with energy 3 and so on. We can easily check itp2d computes these energies correctly.

The simple harmonic potential specified above is also the default potential of itp2d, so simply
running itp2d without any additional parameters will compute the spectrum of the harmonic
oscillator up to the first 16 states. Please note that itp2d will refuse to overwrite an existing
datafile unless you specify the `--force` (or `-f` for short) flag on the command line. Also, to
make itp2d print the energy spectrum after the run is finished you should use the `--verbose` (or
`-v`) flag. So run

	./itp2d -v

(adding the `-f` if necessary) to get the first 16 energies along with the error estimates provided
by the program.

Since the first 16 states are very fast to compute, we can try something larger by increasing the
number of states with the command line flag `-n`. Since higher states of the harmonic oscillator are
spread over a larger area of space, to compute many states we must also make the computation grid
larger. Otherwise the states will become too close to the edge of the simulation box, and they'll
start to notice the periodic boundary conditions. By default itp2d uses a grid of 12 by 12 atomic
units in size, which starts to be too small for 100 states for example. Let's increase the grid size
to 15 units with flag `-l`, and while we're at it, let's also make the grid finer, 100 by 100 grid
points, with the flag `-s`:

	./itp2d -v -n 100 -l 15 -s 100

This still only takes a few seconds. By default itp2d saves it's datafile in `data/itp2d.h5`, you
can change this with the command line flag `-o`. Now we can plot the states we have calculated with
the provided helper script `itp_draw_art.py`. The script requires the [Python Imaging Library][PIL],
so make sure it's installed. Unless you have saved the data elsewhere, simply running

	scripts/itp_draw_art.py

will make a nice coloured 2D plot of the eigenstate densities and save it as a PNG image in
`data/itp2d.png`. Please notice that since the eigenstates of the harmonic oscillator (except the
ground state) are heavily degenerate, you will probably get different eigenstates each time.

[PIL]: http://www.pythonware.com/products/pil

### Solving a particle in a square potential well with a strong magnetic field

To try something more involved, let's look at a particle in a box with a strong external magnetic
field. This is a system with no known analytic solution. We'll start by setting the external
potential to zero and the boundary conditions to Dirichlet with `-p zero --dirichlet`. This will
have the effect of infinite potential walls at the edge of the calculation box. Next we'll set the
box size to pi by pi units with `-l 1.0 --pi`. This is a convenient choice since the spectrum of a
box of this size in *zero* magnetic field is simple. Then we put on a magnetic field with `-B 1.0`.
The value specified is the strength of the magnetic field in SI-bases atomic units. The magnetic field will
be homogeneous and in the direction of the `z`-axis.

When using a magnetic field with Dirichlet boundary conditions there wil be
ringing artifacts near the edges, which arise from the fact that the Dirichlet
boundary conditions are enforced simply by expanding the states in a sine
series instead of a Fourier series. However, the sine functions are
no longer eigenfunctions of the kinetic energy when there is a magnetic field, which causes the
artifacts around the edges of the calculation box. For one, this makes convergence testing based on
the standard deviation of the energy inaccurate, so we need to switch to another test. Let's consider the
states converged in respect to the used time step value when the relative energy change between
successive iterations is less than 0.0001, and completely converged when they converge w.r.t the
timestep with a single iteration. This can be set from the itp2d command line interface with `-T
"relEdelta(0.0001)" -F onestep` (the quotation marks are used only to tell the command shell to not
interpret the parentheses). We'll also need to start with a shorter imaginary time step since the
system is more complicated. To set the time step to 0.1 units we provide `-e 0.1`.

All in all the command line to solve the first 100 states of a particle in a box with a strong
magnetic field is

	./itp2d -f -v -l 1.0 --pi -p zero --dirichlet -B 1.0 -T "relEdelta(0.0001)" -F onestep -n 100 -e 0.1

If you run `itp_draw_art.py` know you'll get some beautiful images. Who knew a particle in a box can
be this beautiful?

### Implementing new potentials

The potential types known to itp2d are defined in the header file `src/potentialtypes.hpp`. You
can create new potential types by creating our own class derived from `PotentialType` and
implementing the required functions. This should be easy to do just by imitating how the existing
potentials are implemented. To include your potential in the in-line documentation you should also edit
`src/commandlineparser.cpp`.

### Reading datafiles created by itp2d with MATLAB

MATLAB has built-in support for reading HDF5 files easily, but the exact syntax
can depend on the version of MATLAB. In version 2009b you can, for example, read the final energies
from a file called `itp2d.h5` with

	hdf5read('itp2d.h5', 'final_energies')

For more in-depth documentation of MATLAB's HDF5 capabilities, please consult
the [MATLAB documentation for the most recent version][matlabhdf5] or [older
versions][matlabhdf5old] (requires login).

[matlabhdf5]: http://www.mathworks.se/help/techdoc/ref/hdf5.html
[matlabhdf5old]: http://www.mathworks.se/help/doc-archives.html

### Accessing itp2d datafiles from Python

There are several libraries for accessing HDF5 data from [Python][] scripts.
The one used in the bundled helper scripts is [h5py][], which imports HDF5 data
as [NumPy][] arrays.

Here is an example of importing a datafile called `itp2d.h5`, listing all the
data it contains, and printing the final energies with h5py:

	import h5py
	from numpy import array
	f = h5py.File("itp2d.h5")
	print f.keys()						# list all datasets in the file
	print array(f["final_energies"])	# read dataset final_energies as an array and print

Please see the helper scripts distributed with itp2d for more examples about
accessing itp2d data from Python.

[h5py]: http://alfven.org/wp/hdf5-for-python/
[NumPy]: http://numpy.scipy.org/

### Providing patches

If you make some changes to itp2d you consider could be of wider use, please send a patch by
uploading you changes to [GitHub][webpage] and sending the author a pull request. You can also simply send
your patch by email directly to the [author][]. All contributions are welcome!
