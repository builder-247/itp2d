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

#include "laplacian.hpp"

// By using FFT to transform into wave-vector space, the Laplacian operator is simply
// a multiplication with -k². The multiplication table of -k² values is saved when the
// operator is created
Laplacian::Laplacian(Transformer const& _tr) : dl(_tr.datalayout), tr(_tr) {
	size_t const& sx = dl.sizex;
	size_t const& sy = dl.sizey;
	double const& normfac = tr.normalization_factor(FFT);
	double kx, ky;
	multipliers = new double[sx*sy];
	for (size_t y=0; y<sy; y++) {
		ky = tr.fft_ky(y);
		for (size_t x=0; x<sx; x++) {
			kx = tr.fft_kx(x);
			dl.value(multipliers, x, y) = (-kx*kx - ky*ky)*normfac;
		}
	}
}

Laplacian::~Laplacian() {
	delete[] multipliers;
}

std::ostream& Laplacian::print(std::ostream& stream) const {
	return stream << "∇²";
}

// Because the multiplication tables were saved in the constructor, the Laplacian is simply
// FFT-multiply-inverseFFT.
void Laplacian::operate(State& state, __attribute__((unused))StateArray& workspace) const {
	assert(dl == state.datalayout);
	state.transform(FFT, tr);
	state.pointwise_multiply(multipliers);
	state.transform(iFFT, tr);
}
