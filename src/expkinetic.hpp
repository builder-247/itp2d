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

#ifndef _EXPKINETIC_HPP_
#define _EXPKINETIC_HPP_

#include "operators.hpp"

// Exponentiated kinetic energy operator for imaginary time propagation.
//
// This code uses the exact factorization of the exponentiated kinetic energy
// operator as specified in M. Aichinger, S. A. Chin and E. Krotscheck,
// Comp. Phys. Comm. 171 (2005), pages 197--207. There is a nice summary of the
// algorithm in page 200, where this kinetic energy part is steps (3)--(5).
// The only difference here is that we use linear gauge instead of symmetric gauge.

class ExpKinetic : public EvolutionOperator {
	public:
		// The constructor creates an operator for
		//
		// 		p·exp(-e·c·T),
		//
		// 	where p is the prefactor, e is the size of the imaginary time step,
		// 	c is some numerical coefficient, and T is the kinetic energy
		// 	operator. Argument B specifies the strength of magnetic field. You
		// 	also need to supply a Transformer instance for doing FFT
		// 	transformations and a boundary type.
		ExpKinetic(double time_step, double B, Transformer const& tr, BoundaryType bt,
				double coeff=-1.0, double prefactor=1.0);
		~ExpKinetic();
		void operate(State& state, __attribute__((unused))StateArray& workspace) const;
		inline size_t required_workspace() const;
		std::ostream& print(std::ostream& out) const;
		inline void set_time_step(double e) { time_step = e; calculate_multipliers(); }
		Transformer const& transformer;
		DataLayout const& datalayout;
		const BoundaryType boundary_type;
		const double B;
		const double coefficient;
		const double prefactor;
	private:
		double time_step;
		// All these multiplier arrays are what the state will be multiplied
		// with after first transforming to a suitable space with FFTs.
		double* multipliers;
		double* xmultipliers;
		double* xmultipliers2;
		double* ymultipliers;
		void calculate_multipliers();
};

inline size_t ExpKinetic::required_workspace() const {
	if (B == 0 or boundary_type == Periodic)
		return 0;
	return 1;
}

#endif // _EXPKINETIC_HPP_
