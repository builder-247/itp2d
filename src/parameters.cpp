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
 * Definitions of defaults
 */

#include "parameters.hpp"

const bool Parameters::default_recover = false;
const unsigned long int Parameters::default_rngseed = 0x20120131;
const char Parameters::default_datafile_name[] = "data/itp2d.h5";
const char Parameters::default_wisdom_file_name[] = "fftw_wisdom";
const Parameters::SaveWhat Parameters::default_save_what = FinalStates;
const bool Parameters::default_clobber = false;
const int Parameters::default_verbosity = 1;
const size_t Parameters::default_num_threads = 2;
const BoundaryType Parameters::default_boundary = Periodic;
const size_t Parameters::default_sizex = 64;
const size_t Parameters::default_sizey = 64;
const double Parameters::default_lenx = 12;
const size_t Parameters::default_N = 25;
const OrthoAlgorithm Parameters::default_ortho_alg = Default;
const Parameters::InitialStatePreset Parameters::default_initialstate_preset = Random;
const char Parameters::default_potential_type[] = "default";
const char Parameters::default_timestep_convergence_test_string[] = "deviation(1e-3,1e-5)";
const char Parameters::default_final_convergence_test_string[] = "deviation(1e-3,0)";
const char Parameters::default_noise_type[] = "none";
const double Parameters::default_B = 0;
const int Parameters::default_halforder = 5;
const double Parameters::default_initial_eps = 0.50;
const double Parameters::default_eps_divisor = 5.0;
const bool Parameters::default_exhaust_eps = false;
const size_t Parameters::default_needed_to_converge = 16;
const size_t Parameters::default_ignore_lowest = 0;
const int Parameters::default_max_steps = 50;
const double Parameters::default_min_time_step = 1e-9;

// Overloaded global operator for easy printing

std::ostream& operator<<(std::ostream& stream, const Parameters& params) {
	stream << "recover: " << params.get_recover() << std::endl;
	stream << "rngseed: " << params.get_random_seed() << std::endl;
	stream << "datafile_name: " << params.get_datafile_name() << std::endl;
	stream << "wisdom_file_name: " << params.get_wisdom_file_name() << std::endl;
	stream << "copy_from: " << params.get_copy_from() << std::endl;
	stream << "save_what: " << params.get_save_what() << std::endl;
	stream << "clobber: " << params.get_clobber() << std::endl;
	stream << "verbosity: " << params.get_verbosity() << std::endl;
	stream << "num_threads: " << params.get_num_threads() << std::endl;
	stream << "ortho_alg: " << params.get_ortho_algorithm() << std::endl;
	stream << "fftw_flags: " << params.get_fftw_flags() << std::endl;
	stream << "sizex: " << params.get_sizex() << std::endl;
	stream << "sizey: " << params.get_sizey() << std::endl;
	stream << "boundary_type: " << params.get_boundary_type() << std::endl;
	stream << "lenx: " << params.get_lenx() << std::endl;
	stream << "N: " << params.get_N() << std::endl;
	stream << "initialstate_description: " << params.get_initialstate_description() << std::endl;
	stream << "potential_type: " << params.get_potential_type() << std::endl;
	stream << "noise: " << params.get_noise().get_description() << std::endl;
	stream << "timestep_convergence_test: " << params.get_timestep_convergence_test().get_description() << std::endl;
	stream << "final_convergence_test: " << params.get_final_convergence_test().get_description() << std::endl;
	stream << "B: " << params.get_B() << std::endl;
	stream << "halforder: " << params.get_halforder() << std::endl;
	stream << "eps_values: " << std::endl;
	std::copy(params.get_eps_values().begin(), params.get_eps_values().end(), std::ostream_iterator<double>(stream, " "));
	stream << "eps_divisor: " << params.get_eps_divisor() << std::endl;
	stream << "exhaust_eps: " << params.get_exhaust_eps() << std::endl;
	stream << "needed_to_converge: " << params.get_needed_to_converge() << std::endl;
	stream << "ignore_lowest: " << params.get_ignore_lowest() << std::endl;
	stream << "max_steps: " << params.get_max_steps() << std::endl;
	stream << "min_time_step: " << params.get_min_time_step() << std::endl;
	return stream;
}

void Parameters::set_to_defaults() {
	// Just set everything to default values
	recover = default_recover;
	rngseed = default_rngseed;
	datafile_name = default_datafile_name;
	wisdom_file_name = default_wisdom_file_name;
	save_what = default_save_what;
	clobber = default_clobber;
	verbosity = default_verbosity;
	num_threads = default_num_threads;
	boundary = default_boundary;
	sizex = default_sizex;
	sizey = default_sizey;
	lenx = default_lenx;
	halforder = default_halforder;
	eps_divisor = default_eps_divisor;
	exhaust_eps = default_exhaust_eps;
	max_steps = default_max_steps;
	min_time_step = default_min_time_step;
	ortho_alg = default_ortho_alg;
	fftw_flags = default_fftw_flags;
	//
	timestep_convergence_test = parse_convergence_description(default_timestep_convergence_test_string);
	final_convergence_test = parse_convergence_description(default_final_convergence_test_string);
	noise = parse_noise_description(default_noise_type);
	define_external_field(default_potential_type, default_B);
	set_num_states(default_N, default_needed_to_converge, default_ignore_lowest);
	define_initial_states(default_initialstate_preset);
}

Parameters::Parameters() :
		timestep_convergence_test(NULL),
		final_convergence_test(NULL) {
	set_to_defaults();
}

Parameters::~Parameters() {
}

void Parameters::define_data_storage(std::string filename, SaveWhat arg_save_what, bool arg_clobber) {
	datafile_name = filename;
	save_what = arg_save_what;
	clobber = arg_clobber;
}

void Parameters::define_grid(size_t arg_sizex, size_t arg_sizey, double arg_lenx, BoundaryType arg_boundary) {
	boundary = arg_boundary;
	sizex = arg_sizex;
	sizey = arg_sizey;
	lenx = arg_lenx;
}

void Parameters::define_external_field(std::string const& ptype, double arg_B) {
	B = arg_B;
	potential_type = ptype;
}

void Parameters::set_num_states(size_t _N, size_t _needed_to_converge, size_t _ignore_lowest) {
	N = _N;
	needed_to_converge = _needed_to_converge;
	ignore_lowest = _ignore_lowest;
	// Sanity check
	if (ignore_lowest + needed_to_converge > N)
		throw InvalidNumberOfStates(N, needed_to_converge, ignore_lowest);
}

void Parameters::define_initial_states(InitialStatePreset preset) {
	initialstate_preset = preset;
	initialstate_func = NULL;
	// Actual implementations are in stateset.hpp because we can't throw around function pointers
	// with references to internal RNG objects.
	switch (preset) {
		case Random:
			initialstate_description = "random";
			break;
		default:
			initialstate_description = "unknown";
			break;
	}
}

void Parameters::define_initial_states(std::string description, initialstatefunc func) {
	initialstate_preset = UserSuppliedInitialState;
	initialstate_func = func;
	initialstate_description = description;
}

void Parameters::set_bailout_limits(int arg_max_steps, double arg_min_time_step) {
	max_steps = arg_max_steps;
	min_time_step = arg_min_time_step;
}
