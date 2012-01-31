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

#ifndef _MULTIPRODUCTSPLIT_HPP_
#define _MULTIPRODUCTSPLIT_HPP_

#include <vector>
#include "operators.hpp"
#include "potential.hpp"
#include "secondordersplit.hpp"
#include "operatorsum.hpp"

// The any-even-order multi-product expansion of the imaginary time evolution operator, as specified in
// S. Chin, arXiv:0809.0914 (2008)

class MultiProductSplit : public EvolutionOperator, public OperatorSum {
	public:
		typedef std::vector<SecondOrderSplit*> membervector;
		MultiProductSplit(int halforder, Potential const& original_potential, double time_step, Transformer const& tr, BoundaryType bt, double B=0);
		~MultiProductSplit();
		void set_time_step(double time_step);
		const int halforder;
	private:
		void calculate_coefficients();
		double* coefficients;
		// The members of the expansion. Each will be a SecondOrderSplit
		// operator with a certain prefactor and time step size
		membervector members;
};

#endif // _MULTIPRODUCTSPLIT_HPP_
