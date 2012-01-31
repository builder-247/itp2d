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

#ifndef _KINETIC_HPP_
#define _KINETIC_HPP_

#include "operators.hpp"

// A kinetic energy operator which possibly includes an external, homogeneous
// magnetic field. The magnetic field is treated in linear gauge, i.e., the
// magnetic vector potential is A = (-By, 0, 0). In all cases the kinetic
// energy operator is applied by using Fourier or sine transforms to expand the
// wave functions in a basis where the kinetic energy is simply a pointwise
// multiplication.

class Kinetic : public virtual Operator {
	public:
		// Construct a kinetic energy operator for given magnetic field strength B.
		Kinetic(double B, Transformer const& tr, BoundaryType bt);
		~Kinetic();
		void operate(State& state, StateArray& workspace) const;
		inline size_t required_workspace() const;
		std::ostream& print(std::ostream& out) const;
		DataLayout const& datalayout;
		Transformer const& transformer;
		const BoundaryType boundary_type;
		const double B;
	private:
		// The multiplication tables
		double* translational_muls_xy;	// If B == 0, the translational part can be done in one go.
		double* translational_muls_x;	// ...but otherwise the multipliers need to be applied
		double* translational_muls_y;	// in two passes.
		double* translational_muls_x2;	// ...and if B != 0 *and* we have Dirichlet boundary we
										// need yet another pass.
		void translation_part(State& state, StateArray& workspace) const;
};

inline size_t Kinetic::required_workspace() const {
	if (B == 0)
		return 0;
	if (boundary_type == Periodic)
		return 1;
	return 2;
}

#endif // _KINETIC_HPP_
