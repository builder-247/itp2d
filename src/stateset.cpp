/* Copyright 2012 Perttu Luukko

 * This file is part of itp2d.

 * itp2d is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.

 * itp2d is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.

 * You should have received a copy of the GNU General Public License along with
 * itp2d.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "stateset.hpp"

// Constructors & Destructors

StateSet::StateSet(size_t arg_N, DataLayout const& dl, OrthoAlgorithm algo) :
		datalayout(dl), N(arg_N), ortho_algorithm(algo), ESolver(N),
		timestep_converged(N), finally_converged(N) {
	dataptr1 = reinterpret_cast<comp*>(fftw_malloc(N*datalayout.N*sizeof(comp)));
	dataptr2 = NULL;
	statearrayptr1 = new StateArray(N, datalayout, dataptr1);
	statearrayptr2 = NULL;
	// More memory is needed if using the HighMem algorithm
	if (ortho_algorithm == HighMem) {
		dataptr2 = reinterpret_cast<comp*>(fftw_malloc(N*datalayout.N*sizeof(comp)));
		statearrayptr2 = new StateArray(N, datalayout, dataptr2);
	}
	state_array = statearrayptr1;
	other_state_array = statearrayptr2;
	overlapmatrix = new comp[N*N];
	for (size_t n=0; n<N; n++) {
		timestep_converged[n] = false;
		finally_converged[n] = false;
	}
	how_many_timestep_converged = 0;
	how_many_finally_converged = 0;
}

StateSet::~StateSet() {
	fftw_free(dataptr1);
	if (dataptr2 != NULL)
		fftw_free(dataptr2);
	delete statearrayptr1;
	delete statearrayptr2;
	delete[] overlapmatrix;
}

// Initializing

void StateSet::init(Parameters const& params, RNG& rng) {
	assert(datalayout.sizex == params.get_sizex());
	assert(datalayout.sizey == params.get_sizey());
	assert(datalayout.dx == params.get_grid_delta());
	Parameters::initialstatefunc func = params.get_initialstate_func();
	// Initialize wave function data
	switch(params.get_initialstate_preset()) {
		case Parameters::UserSuppliedInitialState:
			for (size_t n=0; n<N; n++)
				for (size_t y=0; y<datalayout.sizey; y++) {
					const double dy = datalayout.get_posy(y);
					for (size_t x=0; x<datalayout.sizex; x++) {
						const double dx = datalayout.get_posx(x);
						data(n,x,y) = func(n, dx, dy);
					}
				}
			break;
		case Parameters::CopyFromFile:
			init_from_datafile(params.get_copy_from());
			break;
		case Parameters::Random:
			init_to_gaussian_noise(rng);
			break;
		default:
			throw NotImplemented("Given initial state preset");
			break;
	}
}

void StateSet::init_from_datafile(std::string filename) {
	// open other file read-only
	H5::H5File otherfile;
	otherfile.openFile(filename, H5F_ACC_RDONLY);
	H5::Group otherroot = otherfile.openGroup("/");
	// check that grid properties match
	int othersx, othersy, otherN;
	double otherdx;
	otherroot.openAttribute("num_states").read(H5::PredType::NATIVE_INT, &otherN);
	otherroot.openAttribute("grid_sizex").read(H5::PredType::NATIVE_INT, &othersx);
	otherroot.openAttribute("grid_sizex").read(H5::PredType::NATIVE_INT, &othersy);
	otherroot.openAttribute("grid_delta").read(H5::PredType::NATIVE_DOUBLE, &otherdx);
	if (static_cast<int>(N) != otherN)
		throw GeneralError("Cannot copy state data from datafile: value for num_states does not match.");
	if (static_cast<int>(datalayout.sizex) != othersx)
		throw GeneralError("Cannot copy state data from datafile: value for grid_sizex does not match.");
	if (static_cast<int>(datalayout.sizey) != othersy)
		throw GeneralError("Cannot copy state data from datafile: value for grid_sizey does not match.");
	if (datalayout.dx != otherdx)
		throw GeneralError("Cannot copy state data from datafile: value for grid_delta does not match.");
	// copy data
	H5::DataSet other_states_data = otherfile.openDataSet("/states");
	other_states_data.read(state_array->get_dataptr(), other_states_data.getArrayType());
}

void StateSet::init(comp (*initfunc)(size_t, double, double)) {
	double dx, dy;
	for (size_t n=0; n<N; n++) {
		for (size_t y=0; y<datalayout.sizey; y++) {
			dy = datalayout.get_posy(y);
			for (size_t x=0; x<datalayout.sizex; x++) {
				dx = datalayout.get_posx(x);
				data(n,x,y) = initfunc(n,dx,dy);
			}
		}
	}
}

void StateSet::init_to_gaussian_noise(RNG& rng) {
	for (size_t n=0; n<N; n++)
		for (size_t y=0; y<datalayout.sizey; y++)
			for (size_t x=0; x<datalayout.sizex; x++)
				data(n,x,y) = comp(rng.gaussian_rand(), rng.gaussian_rand());
}

// Orthonormalization with the Löwdin method (also called the subspace orthonormalization method),
// explained for example in M. Aichinger, E. Krotscheck, Comp. Mat. Sci. 34 (2005), pages 193--194.

void StateSet::orthonormalize() throw(std::exception) {
	// Handle the trivial case N=1 separately
	if (N == 1) {
		ortho_timer.start();
		const double norm = (*state_array)[0].norm();
		(*state_array)[0] *= 1.0/norm;
		ortho_timer.stop();
		return;
	}
	ortho_timer.start();
	dot_timer.start();
	// NOTE: Because Eigensolver uses LAPACK, the overlap matrix is stored in column-major format
	#pragma omp parallel for
	for (size_t i=0; i<N; i++) {
		for (size_t j=i; j<N; j++) {
			overlapmatrix[N*j+i] = dot(i,j);
		}
	}
	dot_timer.stop();
	// Solve eigenvalue problem for the overlap matrix
	eigensolve_timer.start();
	ESolver.solve(overlapmatrix);
	for (size_t n=0; n<N; n++) {
		const double eval = ESolver.eigenvalue(n);
		// Check that eigenvalues are OK. If states are propagated "too much",
		// they can become linearly dependent (or close enough so), which
		// causes the overlap matrix to have non-positive eigenvalues and as a
		// result the orthonormalization will fail.
		if (eval <= 0) {
			ortho_timer.stop();
			eigensolve_timer.stop();
			throw(NonPositiveEigenvalue(n, eval, overlapmatrix, N));
		}
		else if (std::fpclassify(eval) != FP_NORMAL) {
			ortho_timer.stop();
			eigensolve_timer.stop();
			throw(NonNormalEigenvalue(n, eval, overlapmatrix, N));
		}
		// Scale eigenvectors with the eigenvalues
		ESolver.scale_eigenvector(overlapmatrix, n, 1/sqrt(eval));
	}
	eigensolve_timer.stop();
	// Form orthonormal states from linear combinations
	lincomb_timer.start();
	const comp one = 1;
	const comp zero = 0;
	const int iN = static_cast<int>(N);
	const int iM = static_cast<int>(datalayout.N);
	comp* const statedata = state_array->get_dataptr();
	switch (ortho_algorithm) {
		case Default:
			// This is the in-place version, which uses less memory
			#pragma omp parallel
			{
				const size_t required_size = N*omp_get_num_threads();
				const size_t thread_offset = N*omp_get_thread_num();
				#pragma omp single
				{
				// Check for enough space on tempstate
				if (tempstate.size() < required_size)
					tempstate.resize(required_size);
				}
				comp* const temp = tempstate.data() + thread_offset;
				comp c;
				#pragma omp for
				for (size_t t=0; t<datalayout.N; t++) {
					// Save old state values
					cblas_zcopy(iN, statedata+t, iM, temp, 1);
					// Note that now overlapmatrix holds the eigenvectors
					cblas_zgemv(CblasRowMajor, CblasNoTrans, iN, iN, &one, overlapmatrix,
							iN, temp, 1, &zero, statedata+t, iM);
				}
			}
			break;
		case HighMem:
			// This is the out-of-place version, where the formation of linear
			// combinations can be expressed simply as a product of two (very
			// large) matrices.
			comp* const other_statedata = other_state_array->get_dataptr();
			assert(statedata != NULL);
			assert(other_statedata != NULL);
			cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, iN, iM, iN, &one,
					overlapmatrix, iN, statedata, iM, &zero, other_statedata, iM);
			switch_state_arrays();
			break;
	}
	lincomb_timer.stop();
	ortho_timer.stop();
}

// Check whether the states are orthonormal to a given precision "epsilon".
bool StateSet::is_orthonormal(double epsilon) const {
	comp z;
	for (size_t n=0; n<N; n++) {
		for (size_t k=0; k<n; k++) {
			z = dot(n,k);
			if (fabs(imag(z)) > epsilon)
				return false;
			if (n != k and fabs(real(z)) > epsilon)
				return false;
			if (n == k and fabs(real(z)-1) > epsilon)
				return false;
		}
	}
	return true;
}

double StateSet::how_orthonormal() const {
	double max = 0;
	for (size_t n=0; n<N; n++) {
		for (size_t k=0; k<=n; k++) {
			const comp z = dot(n,k);
			const double i = imag(z);
			const double r = real(z);
			if (fabs(i) > max)
				max = fabs(i);
			if (n != k) {
				if (fabs(r) > max)
					max = fabs(r);
			}
			else if (fabs(r-1) > max)
				max = fabs(r-1);
		}
	}
	return max;
}
