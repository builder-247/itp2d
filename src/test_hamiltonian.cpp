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

/*
 * Another important test: Compute ground states with ITP and check whether
 * these states are stable in respect to repeated applications of the
 * Hamiltonian operator.
 */

#include "test_hamiltonian.hpp"

class hamiltonian : public testing::Test {
public:
	hamiltonian() : gridsize(40) {
		params.set_random_seed(RNG::produce_random_seed());
		params.set_verbosity(0);
		params.set_fftw_flags(FFTW_ESTIMATE);
		params.define_data_storage("", Parameters::Nothing);
	}
	const size_t gridsize;
	Parameters params;
};

// This is the "skeleton" test function that is applied for different ITPSystems.
void test_hamiltonian_skeleton(ITPSystem* sys, double rmsdist_tolerance,
		__attribute__((unused)) std::string filename) {
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	EXPECT_TRUE(sys->get_error_flag() == false);
	State const& groundstate = sys->get_state(sys->get_sorted_index(0));
	const double groundstateenergy = sys->get_sorted_energy(0);
	OperatorSum const& H = sys->get_hamiltonian();
	Datafile* d = NULL;
	if (dump_data) {
		d = new Datafile(filename, groundstate.datalayout, true);
		d->write_state(0, 0, groundstate);
	}
	// Apply the Hamiltonian several times and save results.
	State temp(groundstate);
	StateArray workspace(3, groundstate.datalayout);
	for (int i=1; i<=3; i++) {
		H(temp, workspace);
		temp /= groundstateenergy;
		// Check that the state has not changed too much.
		const double rmsdist = rms_distance(temp, groundstate);
		EXPECT_LT(rmsdist, rmsdist_tolerance);
		if (dump_data) {
			d->write_state(0, i, temp);
		}
	}
	if (dump_data) {
		delete d;
	}
}

TEST_F(hamiltonian, particle_in_a_box) {
	const double rmsdist_tolerance = 5e-3;
	params.define_grid(gridsize, gridsize, pi, Dirichlet);
	params.define_external_field("zero");
	params.set_num_states(1, 1);
	params.add_eps_value(1.0);
	ITPSystem* sys = new ITPSystem(params);
	test_hamiltonian_skeleton(sys, rmsdist_tolerance, "data/test_hamiltonian_particle_in_a_box.h5");
	delete sys;
}

TEST_F(hamiltonian, harmonic_oscillator) {
	const double rmsdist_tolerance = 5e-3;
	params.define_grid(gridsize, gridsize, 10);
	params.define_external_field("harmonic(1)");
	params.set_num_states(1, 1);
	params.set_operator_splitting_halforder(5);
	params.add_eps_value(1.0);
	ITPSystem* sys = new ITPSystem(params);
	test_hamiltonian_skeleton(sys, rmsdist_tolerance, "data/test_hamiltonian_harmonic.h5");
	delete sys;
}

TEST_F(hamiltonian, harmonic_oscillator_magnetic) {
	const double rmsdist_tolerance = 5e-3;
	params.define_grid(gridsize, gridsize, 10);
	params.define_external_field("harmonic(1)", 1.0);
	params.set_num_states(1, 1);
	params.set_operator_splitting_halforder(5);
	params.add_eps_value(1.0);
	ITPSystem* sys = new ITPSystem(params);
	test_hamiltonian_skeleton(sys, rmsdist_tolerance, "data/test_hamiltonian_harmonic_magnetic.h5");
	delete sys;
}
