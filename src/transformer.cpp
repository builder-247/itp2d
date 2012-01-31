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

#include "transformer.hpp"

// Return the inverse transform corresponding to a transform.
Transform inverse_transform_of(Transform trans) {
	switch (trans) {
		case FFT:
			return iFFT;
		case iFFT:
			return FFT;
		case FFTx:
			return iFFTx;
		case iFFTx:
			return FFTx;
		case FFTy:
			return iFFTy;
		case iFFTy:
			return FFTy;
		case DST:
			return iDST;
		case iDST:
			return DST;
		case DSTx:
			return iDSTx;
		case iDSTx:
			return DSTx;
		case DSTy:
			return iDSTy;
		case iDSTy:
			return DSTy;
		case DCT:
			return iDCT;
		case iDCT:
			return DCT;
		case DCTx:
			return iDCTx;
		case iDCTx:
			return DCTx;
		case DCTy:
			return iDCTy;
		case iDCTy:
			return DCTy;
		default:
			throw GeneralError("Switch statement at inverse_transform_of ended up where it never should.");
			return iFFT;
	}
}

Transformer::Transformer(DataLayout const& lay, unsigned int fftw_flags) :
		datalayout(lay),
		FFT_norm_factor(1.0/static_cast<double>(datalayout.sizex*datalayout.sizey)),
		FFTx_norm_factor(1.0/static_cast<double>(datalayout.sizex)),
		FFTy_norm_factor(1.0/static_cast<double>(datalayout.sizey)),
		DSCT_norm_factor(1.0/static_cast<double>(4*datalayout.sizex*datalayout.sizey)),
		DSCTx_norm_factor(1.0/static_cast<double>(2*datalayout.sizex)),
		DSCTy_norm_factor(1.0/static_cast<double>(2*datalayout.sizey)) {
	const int sx = static_cast<int>(datalayout.sizex);
	const int sy = static_cast<int>(datalayout.sizey);
	const double multiplier_x = M_PI/datalayout.lenx;
	const double multiplier_y = M_PI/datalayout.leny;
	// Compute frequency values
	d_fft_kx = new double[datalayout.sizex];
	d_fft_ky = new double[datalayout.sizey];
	d_dsct_kx = new double[datalayout.sizex];
	d_dsct_ky = new double[datalayout.sizey];
	for (int x=0; x < sx; x++) {
		d_fft_kx[x] = static_cast<double>((x < sx/2)? x : x-sx)*2*multiplier_x;
		d_dsct_kx[x] = static_cast<double>(x+1)*multiplier_x;
	}
	for (int y=0; y < sy; y++) {
		d_fft_ky[y] = static_cast<double>((y < sy/2)? y : y-sy)*2*multiplier_y;
		d_dsct_ky[y] = static_cast<double>(y+1)*multiplier_y;
	}
	// Initialize FFTW plans
	plans = new fftw_plan[num_transform_types];
	// Allocate a temporary data array. This is needed for computing the optimal plans.
	fftw_complex* const fftw_data = reinterpret_cast<fftw_complex*>(fftw_malloc(sx*sy*sizeof(comp)));
	double* const real_data = reinterpret_cast<double*>(fftw_data);
	// plans for plain FFT
	plans[FFT] = fftw_plan_dft_2d(sy, sx, fftw_data, fftw_data, FFTW_FORWARD, fftw_flags);
	plans[iFFT] = fftw_plan_dft_2d(sy, sx, fftw_data, fftw_data, FFTW_BACKWARD, fftw_flags);
	plans[FFTx] = fftw_plan_many_dft(1, &sx, sy, fftw_data, NULL, 1, sx, fftw_data, NULL, 1, sx,
			FFTW_FORWARD, fftw_flags);
	plans[iFFTx] = fftw_plan_many_dft(1, &sx, sy, fftw_data, NULL, 1, sx, fftw_data, NULL, 1, sx,
			FFTW_BACKWARD, fftw_flags);
	plans[FFTy] = fftw_plan_many_dft(1, &sy, sx, fftw_data, NULL, sx, 1, fftw_data, NULL, sx, 1,
			FFTW_FORWARD, fftw_flags);
	plans[iFFTy] = fftw_plan_many_dft(1, &sy, sx, fftw_data, NULL, sx, 1, fftw_data, NULL, sx, 1,
			FFTW_BACKWARD, fftw_flags);
	// For the sine and cosine transforms separately for the real and imaginary
	// part we need to fiddle with the guru interface of FFTW. Some new
	// variables are needed for describing the complex loops involved.
	// Do not try to understand this code without first understanding what fftw_plan_guru does,
	// please refer to the FFTW documentation for that.
	const fftw_r2r_kind DST_kind[] = {FFTW_RODFT10, FFTW_RODFT10};
	const fftw_r2r_kind IDST_kind[] = {FFTW_RODFT01, FFTW_RODFT01};
	const fftw_r2r_kind DCT_kind[] = {FFTW_REDFT10, FFTW_REDFT10};
	const fftw_r2r_kind IDCT_kind[] = {FFTW_REDFT01, FFTW_REDFT01};
	const int n[] = {sy, sx};
	fftw_iodim dimsx[1], dimsy[1], loopsx[2], loopsy[2];
	// x-transform loop setup
	dimsx[0].n = sx;
	dimsx[0].is = 2;
	dimsx[0].os = 2;
	loopsx[0].n = sy;
	loopsx[0].is = 2*sx;
	loopsx[0].os = 2*sx;
	loopsx[1].n = 2;
	loopsx[1].is = 1;
	loopsx[1].os = 1;
	// y-transform loop setup
	dimsy[0].n = sy;
	dimsy[0].is = 2*sx;
	dimsy[0].os = 2*sx;
	loopsy[0].n = sx;
	loopsy[0].is = 2;
	loopsy[0].os = 2;
	loopsy[1].n = 2;
	loopsy[1].is = 1;
	loopsy[1].os = 1;
	// DST plans
	plans[DST] = fftw_plan_many_r2r(2, n, 2, real_data, NULL, 2, 1, real_data, NULL, 2, 1, DST_kind, fftw_flags);
	plans[iDST] = fftw_plan_many_r2r(2, n, 2, real_data, NULL, 2, 1, real_data, NULL, 2, 1, IDST_kind, fftw_flags);
	plans[DSTx] = fftw_plan_guru_r2r(1, dimsx, 2, loopsx, real_data, real_data, DST_kind, fftw_flags);
	plans[iDSTx] = fftw_plan_guru_r2r(1, dimsx, 2, loopsx, real_data, real_data, IDST_kind, fftw_flags);
	plans[DSTy] = fftw_plan_guru_r2r(1, dimsy, 2, loopsy, real_data, real_data, DST_kind, fftw_flags);
	plans[iDSTy] = fftw_plan_guru_r2r(1, dimsy, 2, loopsy, real_data, real_data, IDST_kind, fftw_flags);
	// DCT plans
	plans[DCT] = fftw_plan_many_r2r(2, n, 2, real_data, NULL, 2, 1, real_data, NULL, 2, 1, DCT_kind, fftw_flags);
	plans[iDCT] = fftw_plan_many_r2r(2, n, 2, real_data, NULL, 2, 1, real_data, NULL, 2, 1, IDCT_kind, fftw_flags);
	plans[DCTx] = fftw_plan_guru_r2r(1, dimsx, 2, loopsx, real_data, real_data, DCT_kind, fftw_flags);
	plans[iDCTx] = fftw_plan_guru_r2r(1, dimsx, 2, loopsx, real_data, real_data, IDCT_kind, fftw_flags);
	plans[DCTy] = fftw_plan_guru_r2r(1, dimsy, 2, loopsy, real_data, real_data, DCT_kind, fftw_flags);
	plans[iDCTy] = fftw_plan_guru_r2r(1, dimsy, 2, loopsy, real_data, real_data, IDCT_kind, fftw_flags);
	// free temporary data array
	fftw_free(fftw_data);
}

Transformer::~Transformer() {
	for (size_t i=0; i<num_transform_types; i++)
		fftw_destroy_plan(plans[i]);
	delete[] plans;
	delete[] d_fft_kx;
	delete[] d_fft_ky;
	delete[] d_dsct_kx;
	delete[] d_dsct_ky;
}
