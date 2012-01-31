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

#include "kinetic.hpp"

std::ostream& Kinetic::print(std::ostream& stream) const {
	return stream << "T";
}

Kinetic::Kinetic(double arg_B, Transformer const& tr, BoundaryType bt) :
		datalayout(tr.datalayout), transformer(tr),
		boundary_type(bt), B(arg_B),
		translational_muls_xy(NULL),
		translational_muls_x(NULL),
		translational_muls_y(NULL),
		translational_muls_x2(NULL) {
	size_t const& N = datalayout.N;
	size_t const& sx = datalayout.sizex;
	size_t const& sy = datalayout.sizey;
	double const& normfac = transformer.normalization_factor(boundary_type);
	double const& normfac_x = transformer.normalization_factor_x(boundary_type);
	double const& normfac_y = transformer.normalization_factor_y(boundary_type);
	// Initialize multiplication tables. In all cases the kinetic energy will
	// be calculated by using Fourier transforms to go to a basis where the
	// kinetic energy operator is simply a pointwise multiplication.
	if (B == 0) {
		translational_muls_xy = new double[N];
		for (size_t y=0; y<sy; y++) {
			const double ky = transformer.ky(y, boundary_type);
			for (size_t x=0; x<sx; x++) {
				const double kx = transformer.kx(x, boundary_type);
				datalayout.value(translational_muls_xy, x, y) = 0.5*(kx*kx + ky*ky)*normfac;
			}
		}
	}
	else {
		translational_muls_x = new double[N];
		translational_muls_y = new double[sy];
		switch (boundary_type) {
			case Periodic:
				for (size_t y=0; y<sy; y++) {
					const double ky = transformer.ky(y, boundary_type);
					translational_muls_y[y] = 0.5*ky*ky*normfac_y;
					const double dy = datalayout.get_posy(y);
					for (size_t x=0; x<sx; x++) {
						const double kx = transformer.kx(x, boundary_type) - B*dy;
						datalayout.value(translational_muls_x, x, y) = 0.5*kx*kx*normfac_x;
					}
				}
				break;
			case Dirichlet:
				translational_muls_x2 = new double[N];
				for (size_t y=0; y<sy; y++) {
					const double ky = transformer.ky(y, boundary_type);
					translational_muls_y[y] = 0.5*ky*ky*normfac_y;
					const double dy = datalayout.get_posy(y);
					for (size_t x=0; x<sx; x++) {
						const double kx = transformer.kx(x, boundary_type);
						datalayout.value(translational_muls_x, x, y) = 0.5*(kx*kx+B*B*dy*dy)*normfac_x;
						// In this case the first derivative will give a cosine
						// series, which needs to be multiplied separately
						datalayout.value(translational_muls_x2, x, y) = B*dy*kx*normfac_x;
					}
				}
				break;
		}
	}
}

Kinetic::~Kinetic() {
	delete[] translational_muls_xy;
	delete[] translational_muls_x;
	delete[] translational_muls_y;
	delete[] translational_muls_x2;
}

/*
 * The effect of the kinetic energy operator is calculated by Fourier
 * transforming to wave vector space, where the kinetic energy is simply a
 * multiplication, or a sum of two multiplications in the case of Dirichlet
 * boundary conditions.
 */

void Kinetic::operate(State& state, StateArray& workspace) const {
	if (B == 0) {
		switch (boundary_type) {
			case Periodic:
				state.transform(FFT, transformer);
				state.pointwise_multiply(translational_muls_xy);
				state.transform(iFFT, transformer);
				break;
			case Dirichlet:
				state.transform(DST, transformer);
				state.pointwise_multiply(translational_muls_xy);
				state.transform(iDST, transformer);
				break;
		}
	}
	else {
		State& temp = workspace[0];
		temp = state;
		switch (boundary_type) {
			case Periodic:
				state.transform(FFTx, transformer);
				state.pointwise_multiply(translational_muls_x);
				state.transform(iFFTx, transformer);
				temp.transform(FFTy, transformer);
				temp.pointwise_multiply_y(translational_muls_y);
				temp.transform(iFFTy, transformer);
				break;
			case Dirichlet:
				State& temp2 = workspace[1];
				state.transform(DSTx, transformer);
				temp2 = state;
				state.pointwise_multiply(translational_muls_x);
				state.transform(iDSTx, transformer);
				temp2.pointwise_multiply_imaginary_shiftx(translational_muls_x2);
				temp2.transform(iDCTx, transformer);
				state += temp2;
				temp.transform(DSTy, transformer);
				temp.pointwise_multiply_y(translational_muls_y);
				temp.transform(iDSTy, transformer);
				break;
		}
		state += temp;
	}
}
