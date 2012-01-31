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

#ifndef _LAPLACIAN_HPP_
#define _LAPLACIAN_HPP_

#include "operators.hpp"

// A Laplacian operator. This is used for testing the FFT methods, not really
// for actual ITP simulations

class Laplacian : public virtual Operator {
	public:
		Laplacian(Transformer const& tr);
		~Laplacian();
		void operate(State& state, __attribute__((unused))StateArray& workspace) const;
		inline size_t required_workspace() const { return 0; }
		std::ostream& print(std::ostream& out) const;
	private:
		DataLayout const& dl;
		Transformer const& tr;
		double* multipliers;
};

#endif // _LAPLACIAN_HPP_
