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

#ifndef _POTENTIAL_HPP_
#define _POTENTIAL_HPP_

#include "operators.hpp"
#include "potentialtypes.hpp"
#include "noise.hpp"

// A local potential operator, based on the potential definition given by the PotentialType class, and with
// possibly some added noise.

class Potential : public virtual Operator {
	public:
		Potential(DataLayout const& dl, PotentialType const& ptype, std::string name = "V");
		Potential(DataLayout const& dl, PotentialType const& ptype, Noise const& noise = NoNoise(), std::string name = "V");
		~Potential();
		inline double get_value(size_t x, size_t y) const { return (isnull)? 0 : datalayout.value(values, x, y); }
		inline double const* get_valueptr() const { return (isnull)? NULL : values; }
		inline std::string const& get_name() const { return name; }
		inline void operate(State& state, __attribute__((unused))StateArray& workspace) const;
		inline size_t required_workspace() const { return 0; }
		inline bool is_null() const { return isnull; }
		std::ostream& print(std::ostream& out) const;
		DataLayout const& datalayout;
	private:
		void init_values();
		PotentialType const& type;
		const std::string name;
		bool isnull;
		double* values;
};

// The potential acts simply by multiplying the wave function with the values of the potential, which were
// saved in the contructor
inline void Potential::operate(State& state, __attribute__((unused))StateArray& workspace) const {
	assert(datalayout == state.datalayout);
	if (isnull)
		state.zero();
	else
		state.pointwise_multiply(values);
}

#endif // _POTENTIAL_HPP_
