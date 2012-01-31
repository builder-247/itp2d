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

#include "expkinetic.hpp"

// TODO: A chance for optimization: Switch x<->y in operator factorization so
// that the cheap multiplication is done twice and the expensive one once, not
// the other way around. At this moment this is low-priority, since time spent
// running the propagation gets insignificant when the number of states is
// large.

// For documentation please see the article referenced in expkinetic.hpp.

std::ostream& ExpKinetic::print(std::ostream& stream) const {
	return stream << prefactor << "·exp(" << time_step*coefficient << "·T)";
}

ExpKinetic::ExpKinetic(double e, double B_, Transformer const& tr, BoundaryType bt,
				double c, double p) :
		transformer(tr),
		datalayout(transformer.datalayout),
		boundary_type(bt),
		B(B_),
		coefficient(c),
		prefactor(p),
		time_step(e),
		multipliers(NULL),
		xmultipliers(NULL),
		xmultipliers2(NULL),
		ymultipliers(NULL) {
	if (B == 0) {
		multipliers = new double[datalayout.N];
	}
	else {
		xmultipliers = new double[datalayout.N];
		if (boundary_type == Dirichlet)
			xmultipliers2 = new double[datalayout.N];
		ymultipliers = new double[datalayout.sizey];
	}
	calculate_multipliers();
}

ExpKinetic::~ExpKinetic() {
	delete[] multipliers;
	delete[] xmultipliers;
	delete[] xmultipliers2;
	delete[] ymultipliers;
}

void ExpKinetic::calculate_multipliers() {
	Transformer const& tr = transformer;
	if (B == 0) {
		const double normfac = tr.normalization_factor(boundary_type);
		const double A = coefficient*time_step*0.5;
		const double p = prefactor*normfac;
		for (size_t y=0; y<datalayout.sizey; y++) {
			const double ky = tr.ky(y, boundary_type);
			for (size_t x=0; x<datalayout.sizex; x++) {
				const double kx = tr.kx(x, boundary_type);
				datalayout.value(multipliers, x, y) = p*exp(A*(kx*kx + ky*ky));
			}
		}
	}
	else {
		// Take into account that for periodic boundary conditions the
		// multiplication for the x-part of kinetic energy is done twice, so we
		// need to split the normalization into two
		const double normfacx = (boundary_type == Periodic)?
			sqrt(tr.normalization_factor(FFTx)) : tr.normalization_factor(DSTx);
		const double normfacy = tr.normalization_factor_y(boundary_type);
		const double z = B*time_step;
		// Take care of singularities in the coefficients Cx and Cy for small fields
		const double Cy = (fabs(z) < 1e-6)? 1 + pow(z,2)/6 + pow(z,4)/120			: sinh(z)/z;
		const double Cx = (fabs(z) < 1e-6)? (0.5 + pow(z,2)/24 + pow(z,4)/720)/Cy	: (cosh(z)-1)/(z*sinh(z));
		// Collect constants together
		const double Ax = coefficient*time_step*0.5*Cx;
		const double Ay = coefficient*time_step*0.5*Cy;
		const double px = normfacx;
		const double py = prefactor*normfacy;
		switch (boundary_type) {
			case Periodic:
				for (size_t y=0; y<datalayout.sizey; y++) {
					const double ky = tr.ky(y, boundary_type);
					const double dy = datalayout.get_posy(y);
					ymultipliers[y] = py*exp(Ay*ky*ky);
					for (size_t x=0; x<datalayout.sizex; x++) {
						const double kx = tr.kx(x, boundary_type);
						datalayout.value(xmultipliers, x, y) = px*exp(Ax*(kx-B*dy)*(kx-B*dy));
					}
				}
				break;
			case Dirichlet:
				// For Dirichlet boundaries the situation is more complex:
				// After the exact factorization of the exponentiated kinetic
				// energy operator (see article referenced in expkinetic.hpp),
				// we essentially have to apply an operator of the form
				// exp(c·P²), where c is a constant and P is the kinetic energy
				// operator. When there is a magnetic field, P² will have also
				// first derivatives. Applying this kind of exponentiated first
				// derivative to a sine series will give, after some
				// arithmetic, a series of sines and cosines. This means that
				// the multiplication needs to be split into two, with separate
				// multipliers for the sine and cosine parts of the series.
				for (size_t y=0; y<datalayout.sizey; y++) {
					const double ky = tr.ky(y, boundary_type);
					const double dy = datalayout.get_posy(y);
					const double By = B*dy;
					ymultipliers[y] = py*exp(Ay*ky*ky);
					for (size_t x=0; x<datalayout.sizex; x++) {
						const double kx = tr.kx(x, boundary_type);
						const double commonpart = px*exp(Ax*(kx*kx+By*By));
						const double argument = 2*Ax*kx*By;
						datalayout.value(xmultipliers, x, y) = commonpart*cosh(argument);
						datalayout.value(xmultipliers2, x, y) = commonpart*sinh(argument);
					}
				}
				break;
		}
	}
}

void ExpKinetic::operate(State& state, StateArray& workspace) const {
	assert(datalayout == state.datalayout);
	Transformer const& tr = transformer;
	if (B == 0) {
		// The case with zero magnetic field is simple -- we essentially just
		// multiply with exp(-k²)
		switch (boundary_type) {
			case Periodic:
				state.transform(FFT, tr);
				state.pointwise_multiply(multipliers);
				state.transform(iFFT, tr);
				break;
			case Dirichlet:
				state.transform(DST, tr);
				state.pointwise_multiply(multipliers);
				state.transform(iDST, tr);
				break;
		}
	}
	else {
		switch (boundary_type) {
			case Periodic:
				// For periodic boundaries we can use a Fourier transform to
				// expand the states in plane waves, and the operators will be
				// simple pointwise multiplications in the basis of plane
				// waves. This case is documented well in the article
				// referenced in expkinetic.hpp
				state.transform(FFTx, tr);
				state.pointwise_multiply(xmultipliers);
				state.transform(FFTy, tr);
				state.pointwise_multiply_y(ymultipliers);
				state.transform(iFFTy, tr);
				state.pointwise_multiply(xmultipliers);
				state.transform(iFFTx, tr);
				break;
			case Dirichlet:
				// This is the tricky part: we need to multiply the sine and cosine parts
				// separately
				assert(workspace.size() >= 1);
				State& temp = workspace[0];
				state.transform(DSTx, tr);
				temp = state;
				state.pointwise_multiply(xmultipliers);
				temp.pointwise_multiply_imaginary_shiftx(xmultipliers2);
				state.transform(iDSTx, tr);
				temp.transform(iDCTx, tr);
				state += temp;
				state.transform(DST, tr);
				state.pointwise_multiply_y(ymultipliers);
				state.transform(iDSTy, tr);
				temp = state;
				state.pointwise_multiply(xmultipliers);
				temp.pointwise_multiply_imaginary_shiftx(xmultipliers2);
				state.transform(iDSTx, tr);
				temp.transform(iDCTx, tr);
				state += temp;
				break;
		}
	}
}
