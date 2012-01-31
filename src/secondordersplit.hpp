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

#ifndef _SECONDORDERSPLIT_HPP_
#define _SECONDORDERSPLIT_HPP_

#include "operators.hpp"
#include "potential.hpp"
#include "operatorproduct.hpp"
#include "expkinetic.hpp"
#include "exppotential.hpp"
#include "exceptions.hpp"

// The second order split operator approximation (and powers of it) with or
// without magnetic field. This is the famous Störmer/Verlet factorization
//
//   T_2(e) = exp(-eV/2)·exp(-eT)·exp(-eV/2) + O(e³),
//
// where e is the imaginary time step, V is the (local) potential operator and
// T is the kinetic energy operator.

class SecondOrderSplit : public EvolutionOperator, public OperatorProduct {
	public:
		SecondOrderSplit(Potential const& original_potential, double time_step, double B, Transformer const& tr, BoundaryType bt, double prefactor=1.0, int exponent=1);
		~SecondOrderSplit();
		void set_time_step(double time_step);
	private:
		EvolutionOperator* kinetic_part;
		EvolutionOperator* potential_part;
		EvolutionOperator* potential_part_square;	// Save multiplications by precomputing exp(V)^2
		EvolutionOperator* potential_with_prefactor; // Absorb possible prefactor into this
};

#endif // _SECONDORDERSPLIT_HPP_
