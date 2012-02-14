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

// The parser function; simply uses the generic parse_parameter_string parsers
// and feeds the results to individual subclass constructors
ConvergenceTest const* parse_convergence_description(std::string const& str) {
	// extract name and parameters
	name_parameters_pair p;
	try {
		p = parse_parameter_string(str);
	}
	catch (ParseError& e) {
		std::cerr << e.what() << std::endl;
		throw InvalidPotentialType(str);
	}
	std::string const& name = p.first;
	std::vector<double> const& params = p.second;
	// delegate to constructors based on name
	if (name == "none" or name == "no" or name == "null")
		return new NoConvergenceTest(params);
	if (name == "onestep" or name == "one-step")
		return new OneStepConvergenceTest(params);
	if (name == "absEchange" or name == "absEdelta")
		return new AbsoluteEnergyChangeTest(params);
	if (name == "relEchange" or name == "relEdelta")
		return new RelativeEnergyChangeTest(params);
	if (name == "deviation")
		return new EnergyDeviationChangeTest(params);
	throw UnknownConvergenceType(str);
	return NULL;
}

/* Please see convergence.hpp for documentation of the different test functions */

// NoConvergenceTest

NoConvergenceTest::NoConvergenceTest(std::vector<double> const& params) {
	if (not params.empty())
		throw InvalidConvergenceType("Convergence test NoConvergenceTest does not take parameters");
	init();
}

// OneStepConvergenceTest

OneStepConvergenceTest::OneStepConvergenceTest(std::vector<double> const& params) {
	if (not params.empty())
		throw InvalidConvergenceType("One-step convergence test does not take parameters");
	init();
}

bool OneStepConvergenceTest::test(ITPSystem const& sys, size_t n) const {
	if (sys.get_step_counter() == 1 and sys.get_states().is_timestep_converged(n))
		return true;
	else
		return false;
}

// RelativeEnergyChangeTest

RelativeEnergyChangeTest::RelativeEnergyChangeTest(std::vector<double> const& params) {
	if (params.size() == 1)
		limit = params[0];
	else
		throw InvalidConvergenceType("Convergence test based on relative energy change takes exactly one parameter");
	init();
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

void RelativeEnergyChangeTest::init() {
	std::stringstream ss;
	ss << "relative energy change < " << limit;
	description = ss.str();
}

// AbsoluteEnergyChangeTest

AbsoluteEnergyChangeTest::AbsoluteEnergyChangeTest(std::vector<double> const& params) {
	if (params.size() == 1)
		limit = params[0];
	else
		throw InvalidConvergenceType("Convergence test based on absolute energy change takes exactly one parameter");
	init();
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

void AbsoluteEnergyChangeTest::init() {
	std::stringstream ss;
	ss << "absolute energy change < " << limit;
	description = ss.str();
}

// EnergyDeviationChangeTest

EnergyDeviationChangeTest::EnergyDeviationChangeTest(std::vector<double> const& params) {
	if (params.size() == 1) {
		relative_deviation_limit = params[0];
		difference_limit = 0;
	}
	else if (params.size() == 2) {
		relative_deviation_limit = params[0];
		difference_limit = params[1];
	}
	else
		throw InvalidConvergenceType("Convergence test based on relative standard deviation of energy takes either one or two parameters");
	init();
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

void EnergyDeviationChangeTest::init() {
	std::stringstream ss;
	ss << "relative energy deviation < " << relative_deviation_limit
		<< " or relative energy deviation change < " <<
		difference_limit;
	description = ss.str();
}
