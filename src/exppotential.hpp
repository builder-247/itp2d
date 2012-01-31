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

#ifndef _EXPPOTENTIAL_HPP_
#define _EXPPOTENTIAL_HPP_

#include "operators.hpp"
#include "potential.hpp"

// Exponentiated potential operator. This is simply a local potential operator with exponentiated values.

class ExpPotential : public EvolutionOperator {
	public:
		// This constructor creates a representation for
		//
		//		p·exp(c·e·V),
		//
		//	where p is the prefactor, c is the coefficient, e is imaginary time
		//	step and V is the original potential.
		ExpPotential(Potential const& pot, double time_step, double coefficient=-0.5, double prefactor=1.0);
		~ExpPotential();
		std::ostream& print(std::ostream& out) const;
		inline void operate(State& state, __attribute__((unused))StateArray& workspace) const;
		inline size_t required_workspace() const { return 0; }
		inline void set_time_step(double e) { time_step=e; recalc_potential(); }
		inline void set_coefficient(double coeff) { coefficient=coeff; recalc_potential(); }
		inline void set_prefactor(double prefac) { prefactor=prefac; recalc_potential(); }
		void set_constants(double e, double coeff, double prefac);
		DataLayout const& datalayout;
	private:
		void recalc_potential();
		Potential const& original_potential;
		double* values;
		double time_step;
		double coefficient;
		double prefactor;
		std::string original_name;
		const bool is_trivial; // meaning: is really a unit operator
};

inline void ExpPotential::operate(State& state, __attribute__((unused))StateArray& workspace) const {
	assert(datalayout == state.datalayout);
	if (is_trivial) {
		if (prefactor != 1.0)
			state *= prefactor;
	}
	else {
		state.pointwise_multiply(values);
	}
}

#endif // _EXPPOTENTIAL_HPP_
