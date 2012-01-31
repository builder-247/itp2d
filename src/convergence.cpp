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

#include "convergence.hpp"
// The implementations of convergence tests declared in convergence.hpp require
// the actual class implementation of ITPSystem, so we need to include it here. We
// couldn't include it before since itpsystem.hpp needs to include convergence.hpp.
#include "itpsystem.hpp"

/* Please see convergence.hpp for documentation of the different test functions */

bool OneStepConvergenceTest::test(ITPSystem const& sys, size_t n) const {
	if (sys.get_step_counter() == 1 and sys.get_states().is_timestep_converged(n))
		return true;
	else
		return false;
}

bool RelativeEnergyChangeTest::test(ITPSystem const& sys, size_t n) const {
	if (sys.get_energies().size() < 2)
		return false;
	std::vector<std::vector<double> >::const_iterator it = sys.get_energies().end();
	it--;
	const double thisstep = (*it)[n];
	it--;
	const double prevstep = (*it)[n];
	const double scaled_difference = fabs((thisstep-prevstep)/prevstep);
	return scaled_difference < limit;
}

bool AbsoluteEnergyChangeTest::test(ITPSystem const& sys, size_t n) const {
	if (sys.get_energies().size() < 2)
		return false;
	std::vector<std::vector<double> >::const_iterator it = sys.get_energies().end();
	it--;
	const double thisstep = (*it)[n];
	it--;
	const double prevstep = (*it)[n];
	const double absolute_difference = fabs(thisstep-prevstep);
	return absolute_difference < limit;
}

bool EnergyDeviationChangeTest::test(ITPSystem const& sys, size_t n) const {
	if (sys.get_energies().size() < 2)
		return false;
	std::vector<std::vector<double> > const& energies = sys.get_energies();
	std::vector<std::vector<double> > const& stds = sys.get_standard_deviations();
	const size_t last = stds.size()-1;
	const double thisstep = stds[last][n]/energies[last][n];
	const double prevstep = stds[last-1][n]/energies[last-1][n];
	const bool good = (thisstep < relative_deviation_limit) or (fabs(thisstep - prevstep) < difference_limit);
	return good;
}
