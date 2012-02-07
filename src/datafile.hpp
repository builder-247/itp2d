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
 * A class for storing results of simulations to disk in HDF5 format.
 * The HDF5 C++ interface is mostly C decorated with some classes, but it does the job.
 */

#ifndef _DATAFILE_HPP_
#define _DATAFILE_HPP_

#include <cassert>
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include "H5Cpp.h"
#include "itp2d_common.hpp"
#include "stateset.hpp"
#include "state.hpp"
#include "datalayout.hpp"
#include "potential.hpp"

class Datafile {
	// Some simple constants such as {0, 0} and {1, 1}.
	static const hsize_t unlimited_dims[2];
	static const hsize_t ones[2];
	static const hsize_t zeroes[2];
	public:
		Datafile(std::string filename, DataLayout const& dl, bool clobber = false);
		~Datafile();
		// functions for writing States, StateSets and such into the file
		void write_state(size_t n, size_t m, State const& state);
		void write_stateset(StateSet const& stateset, int step, std::list<size_t> const* sort_order = NULL);
		void write_time_step_history(size_t index, double eps);
		void write_energy_history(std::vector<double> energy_history, size_t index);
		void write_energy_history(std::vector<std::vector<double> > energy_history);
		void write_deviation_history(std::vector<double> energy_history, size_t index);
		void write_deviation_history(std::vector<std::vector<double> > energy_history);
		void write_energies(std::vector<double> energies);
		void write_energy_standard_deviations(std::vector<double> standard_deviations);
		void write_potential(Potential const& pot);
		// functions for adding attributes describing the simulation
		void add_attribute(const char* name, int value);
		void add_attribute(const char* name, unsigned long int value);
		void add_attribute(const char* name, double value);
		void add_attribute(const char* name, const char* value);
		void add_attribute(const char* name, std::string const& value);
		inline void flush() { hfile.flush(H5F_SCOPE_GLOBAL); }
	private:
		static inline void validate_selection(H5::DataSpace const& dataspace);
		void ensure_states_data();
		void ensure_state_history_data();
		void ensure_energies_data();
		void ensure_time_step_history_data();
		void ensure_energy_history_data();
		void ensure_energy_standard_deviations_data();
		void ensure_deviation_history_data();
		void ensure_potential_data();
		DataLayout const& datalayout;
		H5::H5File hfile;
		H5::Group root_group;
		// Datatypes
		H5::DataType const& double_type;
		H5::DataType const& int_type;
		H5::CompType complex_type;
		H5::CompType state_history_type;
		H5::CompType time_step_history_type;
		H5::ArrayType* state_type;		// H5::ArrayType has a protected default
		H5::ArrayType* potential_type;	// constructor, so we need to do this stupid trick.
		// Dataspaces
		const H5::DataSpace scalar_space;
		const H5::DataSpace null_space_1d;
		const H5::DataSpace null_space_2d;
		H5::DataSpace space_1d;
		H5::DataSpace space_2d;
		// Datasets
		H5::DataSet states_data;
		H5::DataSet state_history_data;
		H5::DataSet energies_data;
		H5::DataSet time_step_history_data;
		H5::DataSet energy_history_data;
		H5::DataSet energy_standard_deviations_data;
		H5::DataSet deviation_history_data;
		H5::DataSet potential_data;
		// Dataset property lists
		H5::DSetCreatPropList states_dset_props;
		H5::DSetCreatPropList state_history_props;
		H5::DSetCreatPropList energies_dset_props;
		H5::DSetCreatPropList time_step_history_props;
		H5::DSetCreatPropList energy_history_props;
		H5::DSetCreatPropList potential_dset_props;
};

inline void Datafile::validate_selection(H5::DataSpace const& dataspace) {
	if (!dataspace.selectValid()) {
		hsize_t sel_start[2], sel_end[2];
		dataspace.getSelectBounds(sel_start, sel_end);
		const size_t bbox[4] = {sel_start[0], sel_start[1], sel_end[0], sel_end[1]};
		throw InvalidDataspaceSelection(bbox);
	}
}

// std::pair would be nicer, but unfortunately we need to do this the C way for HDF5.
struct state_history_pair {
	int step;
	int index;
};

struct time_step_history_pair {
	int step;
	double time_step;
};

#endif // _DATAFILE_HPP_
