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
 * The Big Tests: Unit tests for ITPSystem checking few systems with known
 * eigenvalues and eigenstates.
 */

#include "test_itp.hpp"

namespace test_itp_reference {
	// The Fock-Darwin spectrum, i.e., the spectrum of a harmonic oscillator in magnetic field.
	double fock_darwin_energy(int n, int m, double B) {
		const double U = sqrt(1 + 0.25*B*B);
		return (2*n + abs(m) + 1)*U - 0.5*m*B;
	}
	double fock_darwin_state(double r, double B, int n, int m) {
		const double Omega = sqrt(1+0.25*B*B);
		const int absm = abs(m);
		int invC = 1;
		for (int i=1; i<=absm; i++)
			invC *= (n+i);
		const double C = sqrt(2*Omega/(2*pi*static_cast<double>(invC)));
		const double L = std::tr1::assoc_laguerre(n, absm, Omega*r*r);
		return C * exp(-Omega/2*r*r) * pow(Omega,static_cast<double>(absm)/2) * pow(r, absm) * L;
	}
}

class itp : public testing::Test {
public:
	itp() : sx(40), sy(40) {
		params.set_random_seed(RNG::produce_random_seed());
		params.set_verbosity(0);
		params.set_fftw_flags(FFTW_ESTIMATE);
	}
	const size_t sx;
	const size_t sy;
	Parameters params;
};

TEST_F(itp, harmonic_oscillator) {
	const double error_tolerance = 1e-4;
	if (dump_data)
		params.define_data_storage("data/test_itp_harmonic.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(sx, sy, 12.0);
	params.set_num_states(14, 8);
	params.add_eps_value(1.0);
	params.define_external_field("harmonic(1)");
	params.set_final_convergence_test(new EnergyDeviationChangeTest(error_tolerance));
	params.set_timestep_convergence_test(new EnergyDeviationChangeTest(error_tolerance, 0.1*error_tolerance));
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// Compute analytic values of energies
	std::vector<double> reference_energies;
	int E = 1;
	int deg_counter = 1;
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		reference_energies.push_back(E);
		if (deg_counter++ >= E) {
			E++;
			deg_counter = 1;
		}
	}
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		ASSERT_NEAR(sys->get_sorted_energy(n), reference_energies[n], error_tolerance);
	}
	// Cleanup
	delete sys;
}

TEST_F(itp, harmonic_oscillator_dirichlet) {
	const double error_tolerance = 1e-4;
	if (dump_data)
		params.define_data_storage("data/test_itp_harmonic_dirichlet.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(sx, sy, 12.0, Dirichlet);
	params.set_num_states(14, 8);
	params.add_eps_value(1.0);
	params.define_external_field("harmonic(1)");
	params.set_final_convergence_test(new EnergyDeviationChangeTest(error_tolerance));
	params.set_timestep_convergence_test(new EnergyDeviationChangeTest(error_tolerance, 0.1*error_tolerance));
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// Compute analytic values of energies
	std::vector<double> reference_energies;
	int E = 1;
	int deg_counter = 1;
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		reference_energies.push_back(E);
		if (deg_counter++ >= E) {
			E++;
			deg_counter = 1;
		}
	}
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		ASSERT_NEAR(sys->get_sorted_energy(n), reference_energies[n], error_tolerance);
	}
	// Cleanup
	delete sys;
}

TEST_F(itp, harmonic_oscillator_nonsquare_simulation_box) {
	const double error_tolerance = 1e-4;
	if (dump_data)
		params.define_data_storage("data/test_itp_harmonic_nonsquare_simulation_box.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(sx, sy+16, 12.0);
	params.set_num_states(14, 8);
	params.add_eps_value(1.0);
	params.define_external_field("harmonic(1)");
	params.set_final_convergence_test(new EnergyDeviationChangeTest(error_tolerance));
	params.set_timestep_convergence_test(new EnergyDeviationChangeTest(error_tolerance, 0.1*error_tolerance));
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// Compute analytic values of energies
	std::vector<double> reference_energies;
	int E = 1;
	int deg_counter = 1;
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		reference_energies.push_back(E);
		if (deg_counter++ >= E) {
			E++;
			deg_counter = 1;
		}
	}
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		ASSERT_NEAR(sys->get_sorted_energy(n), reference_energies[n], error_tolerance);
	}
	// Cleanup
	delete sys;
}

TEST_F(itp, particle_in_a_box) {
	const double error_tolerance = 1e-4;
	if (dump_data)
		params.define_data_storage("data/test_itp_particle_in_a_box.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(sx, sy, pi, Dirichlet);
	params.set_num_states(16, 12);
	params.define_external_field("zero");
	params.add_eps_value(0.5);
	params.set_final_convergence_test(new EnergyDeviationChangeTest(error_tolerance));
	params.set_timestep_convergence_test(new EnergyDeviationChangeTest(error_tolerance, 0.1*error_tolerance));
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// Compute analytic values of energies
	std::vector<double> reference_energies;
	for (size_t x=1; x<params.get_N(); x++)
		for (size_t y=1; y<params.get_N(); y++)
			reference_energies.push_back(0.5*static_cast<double>(x*x+y*y));
	std::sort(reference_energies.begin(), reference_energies.end());
	reference_energies.resize(params.get_needed_to_converge());
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		ASSERT_NEAR(sys->get_sorted_energy(n), reference_energies[n], error_tolerance);
	}
	// Cleanup
	delete sys;
}

TEST_F(itp, particle_in_a_nonsquare_box) {
	const double error_tolerance = 1e-4;
	if (dump_data)
		params.define_data_storage("data/test_itp_particle_in_a_nonsquare_box.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(2*sx, sy, 2*pi, Dirichlet);
	params.set_num_states(16, 12);
	params.define_external_field("zero");
	params.add_eps_value(0.5);
	params.set_final_convergence_test(new EnergyDeviationChangeTest(error_tolerance));
	params.set_timestep_convergence_test(new EnergyDeviationChangeTest(error_tolerance, 0.1*error_tolerance));
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// Compute analytic values of energies
	std::vector<double> reference_energies;
	for (size_t x=1; x<params.get_N(); x++)
		for (size_t y=1; y<params.get_N(); y++)
			reference_energies.push_back(0.5*(static_cast<double>(x*x)/4+static_cast<double>(y*y)));
	std::sort(reference_energies.begin(), reference_energies.end());
	reference_energies.resize(params.get_needed_to_converge());
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		ASSERT_NEAR(sys->get_sorted_energy(n), reference_energies[n], error_tolerance);
	}
	// Cleanup
	delete sys;
}

TEST_F(itp, particle_in_a_nonsquare_box_otherway) {
	const double error_tolerance = 1e-4;
	if (dump_data)
		params.define_data_storage("data/test_itp_particle_in_a_nonsquare_box_otherway.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(sx, 2*sy, pi, Dirichlet);
	params.set_num_states(16, 12);
	params.define_external_field("zero");
	params.add_eps_value(0.5);
	params.set_final_convergence_test(new EnergyDeviationChangeTest(error_tolerance));
	params.set_timestep_convergence_test(new EnergyDeviationChangeTest(error_tolerance, 0.1*error_tolerance));
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// Compute analytic values of energies
	std::vector<double> reference_energies;
	for (size_t x=1; x<params.get_N(); x++)
		for (size_t y=1; y<params.get_N(); y++)
			reference_energies.push_back(0.5*(static_cast<double>(x*x)/4+static_cast<double>(y*y)));
	std::sort(reference_energies.begin(), reference_energies.end());
	reference_energies.resize(params.get_needed_to_converge());
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		ASSERT_NEAR(sys->get_sorted_energy(n), reference_energies[n], error_tolerance);
	}
	// Cleanup
	delete sys;
}

TEST_F(itp, harmonic_oscillator_magnetic) {
	const double error_tolerance = 1e-3;
	const double B = 10.0;
	if (dump_data)
		params.define_data_storage("data/test_itp_harmonic_magnetic.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(sx, sy, 7.0);
	params.set_num_states(14, 8);
	params.add_eps_value(1.0);
	params.define_external_field("harmonic(1)", B);
	params.set_final_convergence_test(new EnergyDeviationChangeTest(error_tolerance));
	params.set_timestep_convergence_test(new EnergyDeviationChangeTest(error_tolerance, 0.1*error_tolerance));
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// Compute analytic values of energies
	const int N = static_cast<int>(params.get_needed_to_converge());
	std::vector<std::tr1::tuple<double,int,int> > reference_energies;
	for (int n=0; n<N; n++)
		for (int m=-N; m<=N; m++) {
			reference_energies.push_back(std::tr1::make_tuple(test_itp_reference::fock_darwin_energy(n, m, B), n, m));
		}
	std::sort(reference_energies.begin(), reference_energies.end());
	reference_energies.resize(params.get_needed_to_converge());
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		EXPECT_NEAR(sys->get_sorted_energy(n), std::tr1::get<0>(reference_energies[n]), error_tolerance);
	}
	// Compare analytic values of density
	DataLayout const& dl = sys->datalayout;
	Datafile* datafile = NULL;
	State* reference = NULL;
	if (dump_data) {
		reference = new State(dl);
		datafile = new Datafile("data/test_itp_harmonic_magnetic_statecompare.h5", dl, true);
	}
	double sum = 0;
	for (size_t i=0; i<params.get_needed_to_converge(); i++) {
		const int n = std::tr1::get<1>(reference_energies[i]);
		const int m = std::tr1::get<2>(reference_energies[i]);
		State const& S = sys->get_state(sys->get_sorted_index(i));
		for (size_t x=0; x<dl.sizex; x++) {
			const double px = dl.get_posx(x);
			for (size_t y=0; y<dl.sizey; y++) {
				const double py = dl.get_posx(y);
				const double r = hypot(px,py);
				const double ref = test_itp_reference::fock_darwin_state(r, B, n, m);
				if (dump_data)
					(*reference)(x,y) = ref;
				const double refrho = ref*ref;
				const double rho = norm(S(x,y));
				const double rhodiff = rho-refrho;
				sum += rhodiff*rhodiff;
			}
		}
		sum = sqrt(sum/static_cast<double>(dl.N));
		EXPECT_LT(sum, 20*error_tolerance);
		if (dump_data) {
			datafile->write_state(i, 0, S);
			datafile->write_state(i, 1, *reference);
		}
	}
	// Cleanup
	delete sys;
	if (dump_data) {
		delete reference;
		delete datafile;
	}
}

TEST_F(itp, harmonic_oscillator_magnetic_dirichlet) {
	const double error_tolerance = 1e-3;
	const double B = 1.0;
	if (dump_data)
		params.define_data_storage("data/test_itp_harmonic_magnetic_dirichlet.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(sx, sy, 10.0, Dirichlet);
	params.set_num_states(14, 8);
	params.add_eps_value(1.0);
	params.define_external_field("harmonic(1)", B);
	params.set_final_convergence_test(new EnergyDeviationChangeTest(error_tolerance));
	params.set_timestep_convergence_test(new EnergyDeviationChangeTest(error_tolerance, 0.1*error_tolerance));
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// Compute analytic values of energies
	const int N = static_cast<int>(params.get_needed_to_converge());
	std::vector<std::tr1::tuple<double,int,int> > reference_energies;
	for (int n=0; n<N; n++)
		for (int m=-N; m<=N; m++) {
			reference_energies.push_back(std::tr1::make_tuple(test_itp_reference::fock_darwin_energy(n, m, B), n, m));
		}
	std::sort(reference_energies.begin(), reference_energies.end());
	reference_energies.resize(params.get_needed_to_converge());
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		ASSERT_NEAR(sys->get_sorted_energy(n), std::tr1::get<0>(reference_energies[n]), error_tolerance);
	}
	// Compare analytic values of density
	DataLayout const& dl = sys->datalayout;
	Datafile* datafile = NULL;
	State* reference = NULL;
	if (dump_data) {
		reference = new State(dl);
		datafile = new Datafile("data/test_itp_harmonic_magnetic_statecompare.h5", dl, true);
	}
	double sum = 0;
	for (size_t i=0; i<params.get_needed_to_converge(); i++) {
		const int n = std::tr1::get<1>(reference_energies[i]);
		const int m = std::tr1::get<2>(reference_energies[i]);
		State const& S = sys->get_state(sys->get_sorted_index(i));
		for (size_t x=0; x<dl.sizex; x++) {
			const double px = dl.get_posx(x);
			for (size_t y=0; y<dl.sizey; y++) {
				const double py = dl.get_posx(y);
				const double r = hypot(px,py);
				const double ref = test_itp_reference::fock_darwin_state(r, B, n, m);
				if (dump_data)
					(*reference)(x,y) = ref;
				const double refrho = ref*ref;
				const double rho = norm(S(x,y));
				const double rhodiff = rho-refrho;
				sum += rhodiff*rhodiff;
			}
		}
		sum = sqrt(sum/static_cast<double>(dl.N));
		ASSERT_LT(sum, 10*1e-4);
		if (dump_data) {
			datafile->write_state(i, 0, S);
			datafile->write_state(i, 1, *reference);
		}
	}
	// Cleanup
	delete sys;
	if (dump_data) {
		delete reference;
		delete datafile;
	}
}

TEST_F(itp, particle_in_a_box_magnetic) {
	const double error_tolerance = 1e-3;
	if (dump_data)
		params.define_data_storage("data/test_itp_particle_in_a_box_magnetic.h5", Parameters::FinalStates, true);
	else
		params.define_data_storage("", Parameters::Nothing);
	params.define_grid(64, 64, pi, Dirichlet);
	params.set_num_states(16, 8);
	params.define_external_field("zero", 1.0);
	params.add_eps_value(0.3);
	params.set_timestep_convergence_test(new AbsoluteEnergyChangeTest(error_tolerance));
	params.set_final_convergence_test(new OneStepConvergenceTest());
	ITPSystem* sys = new ITPSystem(params);
	while (not sys->is_finished()) {
		sys->step();
	}
	sys->finish();
	ASSERT_FALSE(sys->get_error_flag());
	// No analytic values for energy available, so compare to a "known good" set of values
	double reference_energies[] = {1.077814, 2.134673, 3.089442, 3.551456,
		5.075289, 5.499240, 5.729555, 7.042548, 7.929495, 8.830013, 9.106483,
		9.747584, 10.869846, 11.036034, 12.775542, 13.399877, 13.471929,
		14.310641, 14.751258};
	// Compare
	for (size_t n=0; n<params.get_needed_to_converge(); n++) {
		ASSERT_NEAR(sys->get_sorted_energy(n), reference_energies[n], 2*error_tolerance);
	}
	// Cleanup
	delete sys;
}
