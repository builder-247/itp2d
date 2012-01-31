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

#include "multiproductsplit.hpp"

// For more documentation please see the article referenced in multiproductsplit.hpp.

MultiProductSplit::MultiProductSplit(int hord, Potential const& original_potential, double time_step, Transformer const& tr, BoundaryType bt, double B) :
		halforder((original_potential.is_null())? 1 : hord) {
	assert(tr.datalayout == original_potential.datalayout);
	assert(halforder >= 1);
	coefficients = new double[halforder];
	calculate_coefficients();
	members.resize(halforder);
	for (int t=0; t<halforder; t++) {
		members[t] = new SecondOrderSplit(original_potential, time_step/(t+1), B, tr, bt, coefficients[t], t+1);
		(*this) += *(members[t]);
	}
}

MultiProductSplit::~MultiProductSplit() {
	for (size_t t=0; t<members.size(); t++) {
		delete members[t];
	}
	delete[] coefficients;
}

void MultiProductSplit::set_time_step(double time_step) {
	for (size_t t=0; t<members.size(); t++) {
		members[t]->set_time_step(time_step/static_cast<double>(t+1));
	}
}

void MultiProductSplit::calculate_coefficients() {
	// This is formula (2.12) in the article
	int denom;
	for (int t=1; t<=halforder; t++) {
		denom = 1;
		for (int n=1; n<=halforder; n++) {
			if (n != t)
				denom *= (t*t-n*n);
		}
		coefficients[t-1] = pow(t, 2*(halforder-1))/denom;
	}
}
