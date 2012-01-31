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

/*
 * A class encapsulating discrete {Fourier, sine, cosine} transforms for 2D complex data.
 * Here a sine or cosine transform means the usual real-data transform done *separately* for
 * the real and complex part of the data.
 */

#ifndef _TRANSFORMER_HPP_
#define _TRANSFORMER_HPP_

#include "itp2d_common.hpp"
#include "exceptions.hpp"
#include "datalayout.hpp"

enum Transform { FFT, iFFT, FFTx, iFFTx, FFTy, iFFTy, DST, iDST, DSTx, iDSTx,
	DSTy, iDSTy, DCT, iDCT, DCTx, iDCTx, DCTy, iDCTy };
static const size_t num_transform_types = 18;

Transform inverse_transform_of(Transform trans);

class Transformer {
	public:
		Transformer(DataLayout const& lay, unsigned int fftw_flags = default_fftw_flags);
		~Transformer();
		// operations for querying the frequency values
		inline double const& fft_kx(size_t x) const { return d_fft_kx[x]; }
		inline double const& fft_ky(size_t y) const { return d_fft_ky[y]; }
		inline double const& dsct_kx(size_t x) const { return d_dsct_kx[x]; }
		inline double const& dsct_ky(size_t y) const { return d_dsct_ky[y]; }
		inline double const& kx(size_t x, BoundaryType bt) const
			{ return (bt == Periodic)? d_fft_kx[x] : d_dsct_kx[x]; }
		inline double const& ky(size_t y, BoundaryType bt) const
			{ return (bt == Periodic)? d_fft_ky[y] : d_dsct_ky[y]; }
		inline double const& normalization_factor(Transform trans) const;
		inline double const& normalization_factor(BoundaryType bt) const {
			return (bt==Periodic)? FFT_norm_factor : DSCT_norm_factor;
		}
		inline double const& normalization_factor_x(BoundaryType bt) const {
			return (bt==Periodic)? FFTx_norm_factor : DSCTx_norm_factor;
		}
		inline double const& normalization_factor_y(BoundaryType bt) const {
			return (bt==Periodic)? FFTy_norm_factor : DSCTy_norm_factor;
		}
		// FFT operations
		inline void transform(comp* data, Transform trans) const;
		DataLayout const& datalayout;
	private:
		const double FFT_norm_factor;
		const double FFTx_norm_factor;
		const double FFTy_norm_factor;
		const double DSCT_norm_factor;
		const double DSCTx_norm_factor;
		const double DSCTy_norm_factor;
		double* d_fft_kx;
		double* d_fft_ky;
		double* d_dsct_kx;
		double* d_dsct_ky;
		// FFTW plan structures
		fftw_plan* plans;
};

// Free functions for comparison testing

inline bool operator==(const Transformer& lhs, const Transformer& rhs) {
	if (&lhs == &rhs)
		return true;
	if (lhs.datalayout == rhs.datalayout)
		return true;
	return false;
}

inline bool operator!=(const Transformer& lhs, Transformer const& rhs) {
	return !(lhs == rhs);
}

// Return the normalization factor for the respective transform type.
inline double const& Transformer::normalization_factor(Transform trans) const {
	switch (trans) {
		case FFT:
		case iFFT:
			return FFT_norm_factor;
		case FFTx:
		case iFFTx:
			return FFTx_norm_factor;
		case FFTy:
		case iFFTy:
			return FFTy_norm_factor;
		case DST:
		case iDST:
			return DSCT_norm_factor;
		case DSTx:
		case iDSTx:
			return DSCTx_norm_factor;
		case DSTy:
		case iDSTy:
			return DSCTy_norm_factor;
		case DCT:
		case iDCT:
			return DSCT_norm_factor;
		case DCTx:
		case iDCTx:
			return DSCTx_norm_factor;
		case DCTy:
		case iDCTy:
			return DSCTy_norm_factor;
		default:
			throw GeneralError("Switch statement at Transformer::normalization_factor ended up where it never should.");
			return NaN;
	}
}

// The transform just executes the corresponding FFTW plan. However, for the
// DCT and DST we must first reinterpret the data pointer as a pointer to the
// real part of the first data value.
inline void Transformer::transform(comp* data, Transform trans) const {
	switch (trans) {
		case FFT:
		case iFFT:
		case FFTx:
		case iFFTx:
		case FFTy:
		case iFFTy:
			fftw_execute_dft(plans[trans], data, data);
			break;
		case DST:
		case iDST:
		case DSTx:
		case iDSTx:
		case DSTy:
		case iDSTy:
		case DCT:
		case iDCT:
		case DCTx:
		case iDCTx:
		case DCTy:
		case iDCTy:
			double* rdata = reinterpret_cast<double*>(data);
			fftw_execute_r2r(plans[trans], rdata, rdata);
			break;
	}
}

#endif // _TRANSFORMER_HPP_
