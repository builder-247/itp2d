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

#include "itpsystem.hpp"

// Constructors & destructors

ITPSystem::ITPSystem(Parameters const& given_params,
				volatile sig_atomic_t* arg_abort_flagptr,
				volatile sig_atomic_t* arg_save_flagptr,
				std::ostream& arg_out, std::ostream& arg_err) :
		params(given_params),
		datalayout(params.get_sizex(), params.get_sizey(), params.get_grid_delta()),
		transformer(datalayout, params.get_fftw_flags()),
		boundary_type(params.get_boundary_type()),
		abort_flagptr(arg_abort_flagptr),
		save_flagptr(arg_save_flagptr),
		out(arg_out), err(arg_err),
		finished(false), error_flag(false),
		all_needed_states_timestep_converged(false), all_needed_states_finally_converged(false),
		exhausting_eps_values(params.get_exhaust_eps()),
		rng(params.get_random_seed()),
		pot_type(parse_potential_description(params.get_potential_type())),
		pot(datalayout, *pot_type, rng, params.get_noise()),
		kin(params.get_B(), transformer, boundary_type),
		states(params.get_N(), datalayout, params.get_ortho_algorithm()),
		Esn_tuples(params.get_N()),
		total_step_counter(0),
		step_counter(0), prop_time(0), io_time(0), convtest_time(0),
		eps(NaN), eps_values(params.get_eps_values()) {
	update_timestring();
	omp_set_num_threads(static_cast<int>(params.get_num_threads()));
	if (params.get_save_what() != Parameters::Nothing) {
		// Create a datafile and write some attributes describing the simulation
		datafile = new Datafile(params.get_datafile_name(), datalayout, params.get_clobber());
		datafile->add_attribute("program_version", version_string);
		datafile->add_attribute("random_seed", params.get_random_seed());
		datafile->add_attribute("start_time", timestring);
		datafile->add_attribute("num_threads", params.get_num_threads());
		datafile->add_attribute("num_states", static_cast<int>(params.get_N()));
		datafile->add_attribute("num_wanted_to_converge", static_cast<int>(params.get_needed_to_converge()));
		datafile->add_attribute("ignore_lowest", static_cast<int>(params.get_ignore_lowest()));
		datafile->add_attribute("grid_length", params.get_lenx());
		switch (boundary_type) {
			case Periodic:
				datafile->add_attribute("grid_boundary_type", "periodic");
				break;
			case Dirichlet:
				datafile->add_attribute("grid_boundary_type", "dirichlet");
				break;
		}
		datafile->add_attribute("operator_splitting_order", 2*params.get_halforder());
		datafile->add_attribute("potential", pot_type->get_description());
		datafile->add_attribute("noise", params.get_noise().get_description());
		datafile->add_attribute("timestep_convergence_test", params.get_timestep_convergence_test().get_description());
		datafile->add_attribute("final_convergence_test", params.get_final_convergence_test().get_description());
		datafile->add_attribute("magnetic_field_strength", params.get_B());
		datafile->write_potential(pot);
	}
	else
		datafile = NULL;
	// Start with a imaginary time step (eps) that is either pop'd from a
	// user-supplied list, or if none was supplied, read from the static
	// default definitions of the Parameters class
	if (not eps_values.empty()) {
		eps = eps_values.front();
		eps_values.pop_front();
	}
	else
		eps = Parameters::default_initial_eps;
	if (params.get_save_what() != Parameters::Nothing) {
		datafile->add_attribute("initial_time_step", eps);
		datafile->write_time_step_history(total_step_counter+1, eps);
	}
	// Form the Hamiltonian (sum kinetic and potential energy operators)
	H += kin;
	if (not pot.is_null())
		H += pot;
	// Create an approximation for the imaginary time evolution operator
	T = new MultiProductSplit(params.get_halforder(), pot, eps, transformer, boundary_type, params.get_B());
	// Initialize states
	states.init(params, rng);
	// Allocate some working space for multithreaded operation
	// This is used for operating with the evolution operator
	// and calculating the mean and standard deviation
	const size_t workspace_per_thread = std::max((*T).required_workspace(), H.required_workspace() + 1);
	workslices = new StateArray*[params.get_num_threads()];
	for (size_t i=0; i<params.get_num_threads(); i++) {
		workslices[i] = new StateArray(workspace_per_thread, datalayout);
	}
	if (params.get_save_what() == Parameters::Everything)
		save_states(false);
	io_timer.start();
	if (datafile != NULL)
		datafile->flush();
	io_time += io_timer.stop();
}

ITPSystem::~ITPSystem() {
	delete T;
	for (size_t i=0; i<params.get_num_threads(); i++) {
		delete workslices[i];
	}
	delete[] workslices;
	delete datafile;
	delete pot_type;
}

void ITPSystem::print_initial_message() {
	out << std::fixed << std::setprecision(3) << std::showpoint
		<< "\t"
		<< "converging first " << params.get_needed_to_converge() << " eigenstates starting from " << params.get_ignore_lowest() << ". "
		<< params.get_N() << " states propagated in total" << std::endl
		<< "\tconvergence criteria:" << std::endl
		<< "\t\ttimestep convergence: " << params.get_timestep_convergence_test().get_description() << std::endl
		<< "\t\tfinal convergence: " << params.get_final_convergence_test().get_description() << std::endl
		<< "\tpotential: " << pot_type->get_description() << std::endl
		<< "\t\tnoise: " << params.get_noise().get_description() << std::endl
		<< "\tmagnetic field strength: " << params.get_B() << std::endl
		<< "\tgrid: " << params.get_sizex() << "x" << params.get_sizey() << " of length " << params.get_lenx() << ", ";
	switch (boundary_type) {
		case Periodic:
			out << "periodic boundary conditions" << std::endl;
			break;
		case Dirichlet:
			out << "Dirichlet boundary conditions" << std::endl;
			break;
	}
	if (typeid(*pot_type) == typeid(ZeroPotential)) {
		out << "\tzero potential -> no operator splitting needed" << std::endl;
		assert(T->halforder == 1);
	}
	else
		out << "\toperator splitting order: " << 2*params.get_halforder() << std::endl;
	if (verb(3)) {
		out << "\tsplitted evolution operator: " << *T << std::endl << std::endl
			<< "\traw parameter list: " << std::endl << params << std::endl;
	}
	// Warnings about incompatible or dangerous parameters
	if (boundary_type == Dirichlet and params.get_B() != 0) {
		err << "Warning: you are using Dirichlet boundary conditions with a magnetic field. This can cause slower convergence. Please see the README for more details." << std::endl;
		if (typeid(params.get_timestep_convergence_test()) == typeid(RelativeEnergyDeviationTest) or
			typeid(params.get_final_convergence_test()) == typeid(RelativeEnergyDeviationTest)) {
			err << "Warning: Convergence tests based on the standard deviation of energy become inaccurate when using Dirichlet boundary conditions and a magnetic field. You should use other convergence tests." << std::endl;
		}
	}
}

void ITPSystem::change_time_step() {
	// Pop the next value from the list of time step values, or create new values by dividing with the divisor.
	if (not eps_values.empty()) {
		eps = eps_values.front();
		eps_values.pop_front();
		if (eps_values.empty())
			exhausting_eps_values = false;
	}
	else {
		eps /= params.get_eps_divisor();
	}
	// Bail out if minimum time step is reached
	if (eps < params.get_min_time_step()) {
		err << "Error: Minimum time step reached (" << params.get_min_time_step() << ")." << std::endl
			<< "Bailing out." << std::endl;
		error_flag = true;
		finish();
		return;
	}
	T->set_time_step(eps);
	if (params.get_save_what() != Parameters::Nothing) {
		// record the first step to have the new time step value, i.e. the next value
		// of total_step_counter
		datafile->write_time_step_history(total_step_counter+1, eps);
	}
	step_counter = 0;
	all_needed_states_timestep_converged = false;
	for (size_t n=0; n<params.get_N(); n++)
		states.set_timestep_converged(n, false);
	if (verb(1))
		out << "\tEpsilon changed to " << std::scientific << eps << std::fixed << "." << std::endl;
}

// Check if save_flag is raised by e.g. the signal handlers
void ITPSystem::check_save_flag() {
	if (save_flagptr != NULL and *save_flagptr) {
		save_states();
		*save_flagptr = false;
	}
}

void ITPSystem::propagate() {
	if (verb(2))
		out << "\tPropagating..." << std::endl;
	prop_timer.start();
	#pragma omp parallel for
	for (size_t n=0; n<params.get_N(); n++) {
		// Here we have a chance for optimization, since we could just
		// propagate the non-converged states. However, propagation is a cheap
		// step when the number of states is large, so we'll propagate all
		// states just for added precision and robustness.
		(*T)(states[n], *(workslices[omp_get_thread_num()]));
	}
	prop_time += prop_timer.stop();
}

void ITPSystem::orthonormalize() {
	if (verb(2))
		out << "\tOrthonormalizing..." << std::endl;
	try {
		states.orthonormalize();
	}
	catch (NonPositiveEigenvalue& e) {
		// If the states were propagated too much, i.e., they become too much
		// "almost linearly dependent", orthonormalization fails. If the
		// recovery-mode is on, we can start again by resetting the states and
		// trying a smaller time step value.
		err << "ERROR: " << e.what() << std::endl;
		if (params.get_recover()) {
			err	<< "Trying to recover: Changing time step and resetting states." << std::endl;
			change_time_step();
			states.init(params, rng);
			out << "States reset. Resuming propagation." << std::endl;
		}
		else {
			error_flag = true;
			finish();
		}
		return;
	}
	catch (std::exception& e) {
		err << "ERROR: " << std::endl << e.what() << std::endl;
		err << "Quitting..." << std::endl;
		error_flag = true;
		finish();
		return;
	}
}

void ITPSystem::check_timestep_convergence() {
	convtest_timer.start();
	if (verb(2))
		out << "\tChecking timestep convergence..." << std::endl;
	// Loop through states and check for timestep convergence
	for (size_t n=0; n<params.get_N(); n++) {
		size_t const& index = std::tr1::get<2>(Esn_tuples[n]);
		if (params.get_timestep_convergence_test().test(*this, n))
			states.set_timestep_converged(index, true);
	}
	// Count how many where converged
	bool flag = true;
	size_t count = 0;
	std::vector<size_t> missing;
	for (size_t n = params.get_ignore_lowest(); n<(params.get_ignore_lowest() + params.get_needed_to_converge()); n++) {
		size_t const& index = std::tr1::get<2>(Esn_tuples[n]);
		if (states.is_timestep_converged(index)) {
			count++;
		}
		else {
			flag = false;
			// Keep a record for missing states so that we can tell the user which states we are still waiting for.
			if (missing.size() < 5)
				missing.push_back(n);
		}
	}
	all_needed_states_timestep_converged = flag;
	if (verb(2)) {
		out << "\t\t" << count << "/" << params.get_needed_to_converge()
			<< " wanted states converged in respect to timestep ("
			<< how_many_timestep_converged() << " converged in total)" << std::endl;
		if ((not all_needed_states_timestep_converged)
				and ((params.get_needed_to_converge() - count) <= 5)
				and (not verb(4))) {
			out << "\t\tStill waiting for states: ";
			for (size_t i=0; i<missing.size()-1; i++)
				out << missing[i] << ", ";
			out << missing.back() << std::endl;
		}
	}
	convtest_time += convtest_timer.stop();
}

void ITPSystem::check_final_convergence() {
	convtest_timer.start();
	if (verb(2))
		out << "\tChecking final convergence..." << std::endl;
	for (size_t n=0; n<params.get_N(); n++) {
		size_t const& index = std::tr1::get<2>(Esn_tuples[n]);
		const bool good = params.get_final_convergence_test().test(*this, n);
		states.set_finally_converged(index, good);
	}
	size_t count = 0;
	bool flag = true;
	std::vector<size_t> missing;
	for (size_t n=params.get_ignore_lowest(); n<params.get_ignore_lowest()+params.get_needed_to_converge(); n++) {
		Esn_tuple const& t = Esn_tuples[n];
		const size_t i = std::tr1::get<2>(t);
		if (states.is_finally_converged(i))
			count++;
		else {
			flag = false;
			// Keep a record for missing states so that we can tell the user which states we are still waiting for.
			if (missing.size() < 5)
				missing.push_back(n);
		}
	}
	all_needed_states_finally_converged = flag;
	if (verb(2)) {
		out << "\t\t" << count << "/" << params.get_needed_to_converge() << " wanted states converged ("
			<< how_many_finally_converged() << " converged in total)" << std::endl;
		if ((not all_needed_states_finally_converged)
				and ((params.get_needed_to_converge() - count) <= 5)
				and (not verb(4))) {
			out << "\t\tStill waiting for states: ";
			for (size_t i=0; i<missing.size()-1; i++)
				out << missing[i] << ", ";
			out << missing.back() << std::endl;
		}
	}
	convtest_time += convtest_timer.stop();
}

void ITPSystem::save_states(bool sort) {
	if (verb(2))
		out << "\tSaving states..." << std::endl;
	io_timer.start();
	if (not sort)
		datafile->write_stateset(states, total_step_counter, NULL);
	else {
		// Create a temporary list to store the order of states
		std::list<size_t> sort_order;
		for (size_t n=0; n<params.get_N(); n++)
			sort_order.push_back(std::tr1::get<2>(Esn_tuples[n]));
		datafile->write_stateset(states, total_step_counter, &sort_order);
	}
	io_time += io_timer.stop();
}

void ITPSystem::save_energies() {
	io_timer.start();
	if (not energies.empty())
		datafile->write_energies(energies.back());
	if (not standard_deviations.empty())
		datafile->write_energy_standard_deviations(standard_deviations.back());
	io_time += io_timer.stop();
}

void ITPSystem::save_energy_history() {
	io_timer.start();
	datafile->write_energy_history(energies[total_step_counter-1], total_step_counter-1);
	datafile->write_deviation_history(standard_deviations[total_step_counter-1], total_step_counter-1);
	io_time += io_timer.stop();
}

// A single iteration of imaginary time propagation
void ITPSystem::step() {
	// First check for error conditions and increment some counters
	if (abort_flagptr != NULL and *abort_flagptr) {
		// The abort flag has been raised by a signal handler so we certainly
		// have an error
		error_flag = true;
		finish();
		return;
	}
	if (total_step_counter >= params.get_max_steps()) {
		err << "Error: Maximum number of total steps reached (" << params.get_max_steps() << ")." << std::endl
			<< "Bailing out." << std::endl;
		error_flag = true;
		finish();
		return;
	}
	step_counter++;
	total_step_counter++;
	if (verb(2)) {
		update_timestring();
		out << "Step " << total_step_counter << " (step " << step_counter << " with eps = " << std::scientific << eps << std::fixed << ") starting at "
			<< timestring << std::endl;
	}
	// Propagate the states
	propagate();
	check_save_flag();
	// Orthonormalize them
	orthonormalize();
	check_save_flag();
	// And check convergence
	calculate_energies();
	check_save_flag();
	if (params.get_save_what() != Parameters::Nothing) {
		save_energy_history();
	}
	if (verb(3)) {
		print_energies();
	}
	if (params.get_save_what() == Parameters::Everything) {
		save_states();
	}
	check_timestep_convergence();
	if (all_needed_states_timestep_converged) {
		if (verb(2)) {
			out << "\t\tAll needed states converged in respect to time step with time step = " << std::scientific << eps << std::fixed
				<< " after " << step_counter << " steps." << std::endl;
		}
		check_final_convergence();
		if (all_needed_states_finally_converged) {
			finish();
			return;
		}
		change_time_step();
	}
	else if (exhausting_eps_values)
		change_time_step();
	check_save_flag();
}

// Calculate energies and the standard deviations of energy for each state.
void ITPSystem::calculate_energies() {
	// Don't bother calculating and saving the energies if no steps are done yet.
	if (total_step_counter == 0)
		return;
	if (verb(2))
		out << "\tCalculating energies..." << std::endl;
	convtest_timer.start();
	// Clear the list of (energy,deviation,index)-tuples
	Esn_tuples.clear();
	const size_t N = params.get_N();
	#pragma omp parallel for
	for (size_t n=0; n<N; n++) {
		const std::pair<comp,comp> e_and_sd = H.mean_and_standard_deviation(states[n],
				*(workslices[omp_get_thread_num()]));
		const double energy = std::real(e_and_sd.first);
		const double deviation = std::real(e_and_sd.second);
		const Esn_tuple new_tuple = std::tr1::make_tuple(energy, deviation, n);
		#pragma omp critical
		{
			Esn_tuples.push_back(new_tuple);
		}
	}
	// Sort Esn_tuples according to energy
	sort(Esn_tuples.begin(), Esn_tuples.end());
	// Save energies and standard deviations
	if (energies.size() < static_cast<size_t>(total_step_counter))
		energies.resize(total_step_counter);
	if (standard_deviations.size() < static_cast<size_t>(total_step_counter))
		standard_deviations.resize(total_step_counter);
	std::vector<double>& new_energies = energies[total_step_counter-1];
	std::vector<double>& new_deviations = standard_deviations[total_step_counter-1];
	new_energies.clear();
	new_deviations.clear();
	for (std::vector<Esn_tuple>::const_iterator it = Esn_tuples.begin(); it != Esn_tuples.end(); it++) {
		new_energies.push_back(std::tr1::get<0>(*it));
		new_deviations.push_back(std::tr1::get<1>(*it));
	}
	convtest_time += convtest_timer.stop();
}

void ITPSystem::finish() {
	if (finished)
		return;
	if (params.get_save_what() == Parameters::FinalStates)
		save_states();
	if (params.get_save_what() != Parameters::Nothing) {
		save_energies();
		datafile->add_attribute("num_converged", static_cast<int>(how_many_finally_converged()));
		datafile->add_attribute("error_flag", error_flag);
		datafile->add_attribute("total_steps_done", total_step_counter);
		datafile->add_attribute("propagation_time", get_prop_time());
		datafile->add_attribute("orthonormalization_time", get_ortho_time());
		datafile->add_attribute("dotproduct_time", get_dot_time());
		datafile->add_attribute("eigensolve_time", get_eigensolve_time());
		datafile->add_attribute("lincomb_time", get_lincomb_time());
		datafile->add_attribute("io_time", get_io_time());
		datafile->add_attribute("convtest_time", get_convtest_time());
		datafile->add_attribute("total_time", get_prop_time()+get_ortho_time()+get_convtest_time()+get_io_time());
	}
	finished = true;
	if (verb(2)) {
		update_timestring();
		out << "finish() called at " << timestring << "." << std::endl;
	}
}

void ITPSystem::print_energies() {
	const size_t N = params.get_N();
	out << "\tEnergies:"
		<< std::setprecision(5) << std::showpoint << std::endl;
	for (size_t n=0; n<N; n++) {
		Esn_tuple const& tuple = Esn_tuples[n];
		const double sim_E = std::tr1::get<0>(tuple);
		const double error = std::tr1::get<1>(tuple);
		const size_t index = std::tr1::get<2>(tuple);
		out << "\t" << n << "\t" << std::fixed << sim_E;
		if (not (params.get_B() != 0 and boundary_type == Dirichlet))
			out << " ± " << std::scientific << error;
		if (not states.is_finally_converged(index))
			out << " (not converged)";
		out << std::endl;
	}
}
