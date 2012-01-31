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
 * Unit tests for the Transformer class. This runs discrete Fourier, sine and
 * cosine transforms in cases where the correct result is known and compares
 * the results.
 */

#include "test_transformer.hpp"

namespace test_transformer_reference {
	const int sx = 16;
	const int sy = 16;
	const double dx = 0.5;

	comp initfunc(double x, double y) {
		const double val = exp(-(x*x+y*y));
		return comp(val, -val);
	}

	comp sine_blob(double x, double y) {
		const double val = 4*cos(x*pi/(sx*dx))*cos(y*pi/(sy*dx));
		return comp(val, -val);
	}

	comp unit_wave(double x, double y) {
		const comp i(0,1);
		const double px = 2*pi*x/(sx*dx) - pi/sx;
		const double py = 2*pi*y/(sy*dx) - pi/sy;
		return exp(i*px)*exp(i*py);
	}
}

class transformer : public testing::Test {
public:
	transformer() : dl(test_transformer_reference::sx, test_transformer_reference::sy, test_transformer_reference::dx),
		tr(dl), A(dl, test_transformer_reference::initfunc) {}
	const DataLayout dl;
	const Transformer tr;
	const State A;
};

TEST_F(transformer, kvalues) {
	const int sx = static_cast<int>(dl.sizex);
	const int sy = static_cast<int>(dl.sizey);
	for (int x=0; x<sx; x++) {
		ASSERT_EQ(tr.fft_kx(x), (2*M_PI/dl.lenx)*static_cast<double>((x < sx/2)? x : x-sx));
		ASSERT_EQ(tr.dsct_kx(x), (M_PI/dl.lenx)*static_cast<double>(x+1));
	}
	for (int y=0; y<sy; y++) {
		ASSERT_EQ(tr.fft_ky(y), (2*M_PI/dl.leny)*static_cast<double>((y < sy/2)? y : y-sy));
		ASSERT_EQ(tr.dsct_ky(y), (M_PI/dl.leny)*static_cast<double>(y+1));
	}
}

class transform_type : public testing::TestWithParam<Transform> {
	public:
		transform_type() : dl(test_transformer_reference::sx, test_transformer_reference::sy,
				test_transformer_reference::dx), tr(dl), A(dl, test_transformer_reference::initfunc) {}
		virtual void SetUp() { type = GetParam(); }
		Transform type;
		const DataLayout dl;
		const Transformer tr;
		const State A;
};

TEST_P(transform_type, inverse_transform) {
	State T(A);
	T.transform(type, tr);
	T.transform(inverse_transform_of(type), tr);
	T *= tr.normalization_factor(type);
	ASSERT_LT(rms_distance(A, T), machine_epsilon);
}

INSTANTIATE_TEST_CASE_P(inverse, transform_type, testing::Values(FFT, FFTx, FFTy, DST, DSTx, DSTy, DCT, DCTx, DCTy));

TEST_F(transformer, dirac_delta_idst) {
	State T(dl);
	T.zero();
	T(0,0) = comp(1.0, -1.0);
	T.transform(iDST, tr);
	const State T2(dl, test_transformer_reference::sine_blob);
	ASSERT_LT(rms_distance(T, T2), 10*machine_epsilon);
	if (dump_data) {
		Datafile datafile("data/test_transformer_dirac_delta_idst.h5", dl, true);
		datafile.write_state(0, 0, T);
		datafile.write_state(0, 1, T2);
	}
}

TEST_F(transformer, dirac_delta_ifft) {
	State T(dl);
	T.zero();
	T(1,1) = 1.0;
	T.transform(iFFT, tr);
	const State T2(dl, test_transformer_reference::unit_wave);
	ASSERT_LT(rms_distance(T, T2), 10*machine_epsilon);
	if (dump_data) {
		Datafile datafile("data/test_transformer_dirac_delta_ifft.h5", dl, true);
		datafile.write_state(0, 0, T);
		datafile.write_state(0, 1, T2);
	}
}

TEST_F(transformer, sine_blob_dst) {
	State T(dl, test_transformer_reference::sine_blob);
	T.transform(DST, tr);
	T *= tr.normalization_factor(DST);
	State T2(dl);
	T2.zero();
	T2(0,0) = comp(1.0, -1.0);
	ASSERT_LT(rms_distance(T, T2), 10*machine_epsilon);
	if (dump_data) {
		Datafile datafile("data/test_transformer_sine_blob_dst.h5", dl, true);
		datafile.write_state(0, 0, T);
		datafile.write_state(0, 1, T2);
	}
}

TEST_F(transformer, unit_wave_fft) {
	State T(dl, test_transformer_reference::unit_wave);
	T.transform(FFT, tr);
	T *= tr.normalization_factor(FFT);;
	State T2(dl);
	T2.zero();
	T2(1,1) = 1.0;
	ASSERT_LT(rms_distance(T, T2), 10*machine_epsilon);
	if (dump_data) {
		Datafile datafile("data/test_transformer_unit_wave_fft.h5", dl, true);
		datafile.write_state(0, 0, T);
		datafile.write_state(0, 1, T2);
	}
}
