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
 * A complete simulation program using the imaginary time propagation algorithm
 * to solve the time-independent Schrödinger equation in 2D. This basically
 * just wraps together the ITPSystem class with some command line parameter
 * parsing from CommandLineParser.
 */

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <signal.h>

#include "parameters.hpp"
#include "commandlineparser.hpp"
#include "itpsystem.hpp"
#include "timer.hpp"

using namespace std;

// We need to have a global pointer to call methods from the signal handlers.
// This is ugly and dangerous, but I see no way around it.
ITPSystem* sys = NULL;

// Signal handler for SIGINT and SIGUSR1
void sighandler(int s) {
	if (s == SIGINT) {
		if (sys == NULL or sys->get_quit_flag()) {
			exit(1);
		}
		cerr << "Caught SIGINT. Saving data and quitting at next convenient spot." << endl
			<< "Press Ctrl-C again to signal immediate stop." << endl;
		sys->emergency_quit();
	}
	if (s == SIGUSR1) {
		cerr << "Caught SIGUSR1. Saving states at next convenient stop and continuing." << endl;
		sys->save_states_at_next_opportunity();
	}
}

int main(int argc, char* argv[]) {
	Timer timer;
	RNG rng;
	// Trap SIGINT and SIGUSR1
	signal(SIGINT, sighandler);
	signal(SIGUSR1, sighandler);
	// Parse parameters
	vector<string> args(argv, argv+argc);
	CommandLineParser parser;
	parser.parse(args);
	Parameters params(parser.get_params());
	params.define_random_source(rng);
	// Import FFTW Wisdom if available
	std::string const& fftw_wisdom_filename = params.get_wisdom_file_name();
	FILE* wisdom_file = fopen(fftw_wisdom_filename.c_str(), "r");
	if (wisdom_file != NULL) {
		fftw_import_wisdom_from_file(wisdom_file);
		fclose(wisdom_file);
	}
	// Initialize ITPSystem
	if (params.get_verbosity() >= 1) {
		cout << "Initializing ITP system..." << endl;
	}
	sys = new ITPSystem(params);
	if (params.get_verbosity() >= 1) {
		sys->print_initial_message();
		cout << "Initializations ready. Starting propagation." << endl;
		if (params.get_verbosity() >= 2)
			cout << endl;
	}
	timer.start();
	/*** Main loop starts ***/
	while(not sys->is_finished()) {
		sys->step();
	}
	/*** Main loop ends ***/
	const double elapsed_time = timer.stop();
	const bool error_flag = sys->get_error_flag();
	if (params.get_verbosity() >= 1) {
		if (params.get_verbosity() >= 2)
			cout << endl;
		cout << "Finished." << endl;
		cout << "Total " << sys->get_total_step_counter() << " steps of ITP performed." << endl;
		cout << fixed << "Total simulation time: " << elapsed_time << " s." << endl;
	}
	if (params.get_verbosity() >= 2) {
		const double prop_ratio = sys->get_prop_time()/elapsed_time;
		const double ortho_ratio = sys->get_ortho_time()/elapsed_time;
		const double dot_ratio = sys->get_dot_time()/elapsed_time;
		const double eigensolve_ratio = sys->get_eigensolve_time()/elapsed_time;
		const double lincomb_ratio = sys->get_lincomb_time()/elapsed_time;
		const double convtest_ratio = sys->get_convtest_time()/elapsed_time;
		const double io_ratio = sys->get_io_time()/elapsed_time;
		const double other_ratio = 1.0 - prop_ratio - ortho_ratio - convtest_ratio - io_ratio;
		cout << "Ratios:" << std::fixed << std::setprecision(3) << endl
			<< "\tPropagation:         " << prop_ratio << endl
			<< "\tOrthonormalization:  " << ortho_ratio << endl
			<< "\t      dot products:  " << dot_ratio << endl
			<< "\t      eigenvalues:   " << eigensolve_ratio << endl
			<< "\t      combination:   " << lincomb_ratio << endl
			<< "\tConvergence testing: " << convtest_ratio << endl
			<< "\tI/O:                 " << io_ratio << endl
			<< "\tOther:               " << other_ratio << endl;
		cout << endl;
		if (not sys->get_error_flag())
			sys->print_energies();
	}
	// Save FFTW Wisdom
	wisdom_file = fopen(fftw_wisdom_filename.c_str(), "w");
	fftw_export_wisdom_to_file(wisdom_file);
	delete sys;
	fftw_cleanup();
	if (error_flag)
		return 1;
	else
		return 0;
}
