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
 * A class for storing data on HDF5 files easily.
 */

#include "datafile.hpp"

// Defining trivial constants seems very stupid indeed but we need to do this
// since many HDF5 methods require pointers to arrays, and we can't have
// pointers to anonymous constants.
const hsize_t Datafile::unlimited_dims[] = {H5S_UNLIMITED, H5S_UNLIMITED};
const hsize_t Datafile::ones[] = {1, 1};
const hsize_t Datafile::zeroes[] = {0, 0};

Datafile::Datafile(std::string filename, DataLayout const& dl, bool clobber) :
		datalayout(dl),
		double_type(H5::PredType::NATIVE_DOUBLE),
		int_type(H5::PredType::NATIVE_INT),
		scalar_space(H5S_SCALAR),
		null_space_1d(1, zeroes, unlimited_dims),
		null_space_2d(2, zeroes, unlimited_dims) {
	H5::Exception::dontPrint(); // Turn off error printing from HDF5. Use exceptions.
	const hsize_t state_dims[2] = {datalayout.sizex, datalayout.sizey};
	// Try to open the file, overwriting it if clobber is set
	try {
		hfile = H5::H5File(filename, clobber? H5F_ACC_TRUNC : H5F_ACC_EXCL);
	}
	catch (H5::Exception& e) {
		std::cerr << "Caught HDF5 error when trying to create datafile '" << filename << "'" << std::endl
			<< "Did you try writing on an existing datafile, or write to" << std::endl
			<< "a non-existant directory?" << std::endl
			<< "Detailed error report from HDF5:" << std::endl << std::endl;
		e.printError();
		throw;
	}
	root_group = hfile.openGroup("/");
	// Datatypes
	complex_type = H5::CompType(sizeof(comp));
	complex_type.insertMember("r", 0, double_type);
	complex_type.insertMember("i", sizeof(double), double_type);
	// if states are written to the file several times, state_history will
	// record the iteration (or step) the states were written on, and the index
	// the states were saved to
	state_history_type = H5::CompType(2*sizeof(int));
	state_history_type.insertMember("step", offsetof(state_history_pair, step), int_type);
	state_history_type.insertMember("index", offsetof(state_history_pair, index), int_type);
	// time_step_history does the same, but it saves the current time step size
	time_step_history_type = H5::CompType(8 + sizeof(double));
	time_step_history_type.insertMember("step", offsetof(time_step_history_pair, step), int_type);
	time_step_history_type.insertMember("time_step", offsetof(time_step_history_pair, time_step), double_type);
	state_type = new H5::ArrayType(complex_type, 2, state_dims);
	potential_type = new H5::ArrayType(double_type, 2, state_dims);
	// Dataspaces
	// Everything is just initialized to zero size and expanded from there
	space_1d = H5::DataSpace(1, zeroes, unlimited_dims);
	space_2d = H5::DataSpace(2, zeroes, unlimited_dims);
	// Dataset properties, these turn on inline compression of data
	states_dset_props.setChunk(2, ones);
	states_dset_props.setDeflate(9);
	state_history_props.setChunk(1, ones);
	state_history_props.setDeflate(9);
	energies_dset_props.setChunk(1, ones);
	energies_dset_props.setDeflate(9);
	time_step_history_props.setChunk(1, ones);
	time_step_history_props.setDeflate(9);
	energy_history_props.setChunk(2, ones);
	energy_history_props.setDeflate(9);
	// Standard attributes
	add_attribute("grid_sizex", static_cast<int>(datalayout.sizex));
	add_attribute("grid_sizey", static_cast<int>(datalayout.sizey));
	add_attribute("grid_delta", datalayout.dx);
}

Datafile::~Datafile() {
	delete state_type;
	delete potential_type;
	hfile.close();
}

/*
 * The ensure-functions just check if the dataset exists and create it if it
 * doesn't. This lazy initialization trick is useful, since not all simulations
 * save all possible data.
 */

void Datafile::ensure_states_data() {
	if (states_data.getId() == 0)
		states_data = hfile.createDataSet("/states",
				*state_type, null_space_2d, states_dset_props);
}

void Datafile::ensure_state_history_data() {
	if (state_history_data.getId() == 0)
		state_history_data = hfile.createDataSet("/state_history",
				state_history_type, null_space_1d, state_history_props);
}

void Datafile::ensure_energies_data() {
	if (energies_data.getId() == 0)
		energies_data = hfile.createDataSet("/final_energies",
				double_type, null_space_1d, energies_dset_props);
}

void Datafile::ensure_time_step_history_data() {
	if (time_step_history_data.getId() == 0)
		time_step_history_data = hfile.createDataSet("/time_step_history",
				time_step_history_type, null_space_1d, time_step_history_props);
}

void Datafile::ensure_energy_history_data() {
	if (energy_history_data.getId() == 0)
		energy_history_data = hfile.createDataSet("/energy_history",
				double_type, null_space_2d, energy_history_props);
}

void Datafile::ensure_energy_standard_deviations_data() {
	if (energy_standard_deviations_data.getId() == 0)
		energy_standard_deviations_data = hfile.createDataSet("/final_energy_standard_deviations",
				double_type, null_space_1d, energies_dset_props);
}

void Datafile::ensure_deviation_history_data() {
	if (deviation_history_data.getId() == 0)
		deviation_history_data = hfile.createDataSet("/deviation_history",
				double_type, null_space_2d, energy_history_props);
}

void Datafile::ensure_potential_data() {
	if (potential_data.getId() == 0) {
		potential_data = hfile.createDataSet("/potential_values", *potential_type, scalar_space, potential_dset_props);
	}
}

// Functions for writing data to the file

// States are written in a 2D array. The second index m is meant to refer to
// the index of the state in the set of states, whereas the first index is a
// running index for saving the same set of states several times, e.g., after
// each iteration
void Datafile::write_state(size_t n, size_t m, State const& state) {
	assert(datalayout == state.datalayout);
	const hsize_t coords[1][2] = {{n, m}};
	const hsize_t min_size[2] = {n+1, m+1};
	try {
		// extend dataset if needed, and write the state at the position specified by indices n and m
		ensure_states_data();
		states_data.extend(min_size);
		space_2d.setExtentSimple(2, min_size);
		space_2d.selectElements(H5S_SELECT_SET, 1, reinterpret_cast<const hsize_t*>(coords));
		validate_selection(space_2d);
		states_data.write(state.data_ptr(), *state_type, scalar_space, space_2d);
	}
	catch(H5::Exception& e) {
		e.printError();
		throw;
	}
}


// Write a whole StateSet in one go. The optional argument sort_order can be
// given to write states in an order different from their order in the
// StateSet (for example, in order of increasing energy).
void Datafile::write_stateset(StateSet const& stateset, int step, std::list<size_t> const* sort_order) {
	hsize_t states_cur_size[2];
	hsize_t state_history_cur_size[1];
	try {
		ensure_states_data();
		ensure_state_history_data();
		// calculate current size of states_data
		states_data.getSpace().getSimpleExtentDims(states_cur_size);
		// new slot for states will be same as states_cur_size[0]
		const hsize_t new_slot = states_cur_size[0];
		state_history_pair pair;
		pair.step = step;
		pair.index = static_cast<int>(new_slot);
		// calculate current size of state_history_data and extend by one
		state_history_data.getSpace().getSimpleExtentDims(state_history_cur_size);
		hsize_t new_size = state_history_cur_size[0]+1;
		state_history_data.extend(&new_size);
		space_1d.setExtentSimple(1, &new_size);
		space_1d.selectElements(H5S_SELECT_SET, 1, state_history_cur_size);
		validate_selection(space_1d);
		// record what was the step when states were saved
		state_history_data.write(&pair, state_history_type, scalar_space, space_1d);
		// write states
		const size_t N = stateset.get_num_states();
		if (sort_order == NULL)
			for (size_t m=0; m<N; m++)
				write_state(new_slot, m, stateset[m]);
		else {
			std::list<size_t>::const_iterator it = sort_order->begin();
			for (size_t m=0; m<N; m++) {
				write_state(new_slot, m, stateset[*it]);
				it++;
			}
		}
	}
	catch(H5::Exception& e) {
		e.printError();
		throw;
	}
}

void Datafile::write_time_step_history(size_t index, double eps) {
	// pack values of index and eps into a struct in memory
	time_step_history_pair pair;
	pair.step = static_cast<int>(index);
	pair.time_step = eps;
	try {
		ensure_time_step_history_data();
		// calculate current size, extend if needed, and write data. You know the drill.
		H5::DataSpace tempspace = time_step_history_data.getSpace();
		tempspace.selectAll();
		hsize_t cur_size = tempspace.getSelectNpoints();
		hsize_t new_size = cur_size + 1;
		space_1d.setExtentSimple(1, &new_size);
		time_step_history_data.extend(&new_size);
		space_1d.selectElements(H5S_SELECT_SET, 1, &cur_size);
		validate_selection(space_1d);
		time_step_history_data.write(&pair, time_step_history_type, scalar_space, space_1d);
	}
	catch(H5::Exception& e) {
		e.printError();
		throw;
	}
}

void Datafile::write_energy_history(std::vector<double> energy_history, size_t index) {
	const hsize_t N = energy_history.size();
	const hsize_t min_size[2] = {index+1, N};
	const hsize_t count[2] = {1, N};
	const hsize_t start[2] = {index, 0};
	try {
		ensure_energy_history_data();
		energy_history_data.extend(min_size);
		space_2d.setExtentSimple(2, min_size);
		space_2d.selectHyperslab(H5S_SELECT_SET, count, start);
		validate_selection(space_2d);
		space_1d.setExtentSimple(1, &N);
		space_1d.selectAll();
		energy_history_data.write(&(energy_history.front()), double_type,
				space_1d, space_2d);
	}
	catch (H5::Exception& e) {
		e.printError();
		throw;
	}
}

void Datafile::write_energy_history(std::vector<std::vector<double> > energy_history) {
	const hsize_t N = energy_history.size();
	for (size_t n=0; n<N; n++)
		write_energy_history(energy_history[n], n);
}

void Datafile::write_energies(std::vector<double> energies) {
	const hsize_t N = energies.size();
	try {
		ensure_energies_data();
		space_1d.setExtentSimple(1, &N);
		space_1d.selectAll();
		energies_data.extend(&N);
		energies_data.write(&energies.front(), double_type, space_1d, space_1d);
	}
	catch (H5::Exception& e) {
		e.printError();
		throw;
	}
}

void Datafile::write_deviation_history(std::vector<double> deviation_history, size_t index) {
	const hsize_t N = deviation_history.size();
	const hsize_t min_size[2] = {index+1, N};
	const hsize_t count[2] = {1, N};
	const hsize_t start[2] = {index, 0};
	try {
		ensure_deviation_history_data();
		deviation_history_data.extend(min_size);
		space_2d.setExtentSimple(2, min_size);
		space_2d.selectHyperslab(H5S_SELECT_SET, count, start);
		validate_selection(space_2d);
		space_1d.setExtentSimple(1, &N);
		space_1d.selectAll();
		deviation_history_data.write(&(deviation_history.front()), double_type,
				space_1d, space_2d);
	}
	catch (H5::Exception& e) {
		e.printError();
		throw;
	}
}

void Datafile::write_deviation_history(std::vector<std::vector<double> > deviation_history) {
	const hsize_t N = deviation_history.size();
	for (size_t n=0; n<N; n++)
		write_deviation_history(deviation_history[n], n);
}

void Datafile::write_energy_standard_deviations(std::vector<double> standard_deviations) {
	const hsize_t N = standard_deviations.size();
	try {
		ensure_energy_standard_deviations_data();
		space_1d.setExtentSimple(1, &N);
		space_1d.selectAll();
		energy_standard_deviations_data.extend(&N);
		energy_standard_deviations_data.write(&standard_deviations.front(), double_type, space_1d, space_1d);
	}
	catch (H5::Exception& e) {
		e.printError();
		throw;
	}
}

void Datafile::write_potential(Potential const& pot) {
	if (pot.is_null())
		return;	// Potential is zero so we have nothing to write.
	try {
		ensure_potential_data();
		potential_data.write(pot.get_valueptr(), *potential_type, scalar_space, scalar_space);
	}
	catch(H5::Exception& e) {
		e.printError();
		throw;
	}
}

/*
 * The attibute adding functions are just decorated overloaded wrappers over createAttribute.
 */

void Datafile::add_attribute(const char* name, double value) {
	try {
		H5::Attribute attr = root_group.createAttribute(name, double_type, scalar_space);
		attr.write(double_type, &value);
	}
	catch (H5::Exception& e) {
		e.printError();
		throw;
	}
}

void Datafile::add_attribute(const char* name, int value) {
	try {
		H5::Attribute attr = root_group.createAttribute(name, int_type, scalar_space);
		attr.write(int_type, &value);
	}
	catch (H5::Exception& e) {
		e.printError();
		throw;
	}
}

void Datafile::add_attribute(const char* name, const char* value) {
	add_attribute(name, std::string(value));
}

void Datafile::add_attribute(const char* name, std::string const& value) {
	const size_t len = value.length();
	try {
		H5::AtomType string_type = H5::AtomType(H5::PredType::C_S1);
		// setSize does not like zero sizes
		if (len != 0)
			string_type.setSize(len);
		H5::Attribute attr = root_group.createAttribute(name, string_type, scalar_space);
		attr.write(string_type, value.c_str());
	}
	catch (H5::Exception& e) {
		e.printError();
		throw;
	}
}
