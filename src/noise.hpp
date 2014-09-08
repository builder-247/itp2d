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

#ifndef _NOISE_HPP_
#define _NOISE_HPP_

#include <list>
#include <vector>
#include <tr1/tuple>

#include "itp2d_common.hpp"
#include "datalayout.hpp"
#include "rng.hpp"
#include "parser.hpp"
#include "constraint.hpp"

// Public noise class interface

class Noise {
	public:
		virtual ~Noise() {};
		virtual void add_noise(DataLayout const& dl, double* pot_values) const = 0;
		virtual std::string const& get_description() const = 0;
		// Write internal data which can be used to re-create the noise realization. For example for
		// Gaussian impurities, store the positions, amplitudes and widths of the Gaussians.
		virtual void write_realization_data(std::vector<double>& vec) const = 0;
};

// Interface class of impurity types

class ImpurityType {
	public:
		virtual ~Impurity() {};
		virtual void new_realization(std::vector<double>& params) = 0;
		virtual void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const = 0;
		std::string const& get_description() const { return description; }
	private:
		std::string description;
};

// Interface class of impurity distributions

class ImpurityDistribution {
	public:
		typedef std::pair<double, double> coordinate_pair;
		virtual ~ImpurityDistribution() {};
		std::list<coordinate_pair> const& get_coordinates() const { return coordinates; }
		std::string const& get_description() const { return description; }
	private:
		std::string description;
		std::list<coordinate_pair> coordinates;
};

// Regular spatial noise class built from impurity type + impurity distribution

class SpatialImpurities : public Noise {
	public:
		SpatialImpurities(ImpurityType const& type, ImpurityDistribution& distribution);
		virtual void add_noise(DataLayout const& dl, double* pot_values) const;
		std::string const& get_description() const;
		void write_realization_data(std::vector<double>& vec) const;
	private:
		ImpurityType const& type;
		ImpurityDistribution const& distribution;
		std::list<double> realization_data;
};

// The trivial case of no noise at all

class NoNoise : public Noise {
	public:
		void add_noise(__attribute__((unused)) DataLayout const& dl, __attribute__((unused)) double* pot_values) const {}
		static std::string const& get_description() const { return "none"; }
		void write_realization_data(std::vector<double>& vec) const { vec.clear(); }
};

// Individual impurity types

class GaussianImpurities : public ImpurityType {
	public:
		GaussianImpurities(double amp_mean, double amp_stdev, double width_mean, double width_stdev, RNG& rng);
		void new_realization(std::vector<double>& params);
		void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const;
		std::string const& get_description() const;
	private:
		double amp_mean;
		double amp_stdev;
		double width_mean;
		double width_stdev;
		RNG& rng;
};

// Individual distribution types


// A parser function for returning a Noise instance from a user provided
// desciption string

Noise const* parse_noise_description(std::string const& str, DataLayout const& dl, Constraint const& constraint, RNG& rng);

// Individual noise types

/* The simple case of no noise at all */

/* Randomly distributed gaussian spikes with normally distributed width
 * and normally distributed amplitude. */
class GaussianNoise : public Noise {
	public:
		typedef std::tr1::tuple<double, double, double, double> spike; // x-coord, y-coord, amplitude and width of a spike
		GaussianNoise(double d, double amp, double amp_stdev, double width, double width_stdev,
				DataLayout const& dl, Constraint const& constr, RNG& rng);
		void add_noise(DataLayout const& dl, double* pot_values) const;
		void write_realization_data(std::vector<double>& vec) const;
		DataLayout const& datalayout;
		Constraint const& constraint;
	private:
		double density;
		double amplitude_mean;
		double amplitude_stdev;
		double width_mean;
		double width_stdev;
		void init(RNG& rng);
		std::vector<spike> spikes;
};

class SingleGaussianNoise : public Noise {
	public:
		SingleGaussianNoise(double x, double y, double amplitude, double width);
		void add_noise(DataLayout const& dl, double* pot_values) const;
		void write_realization_data(std::vector<double>& vec) const;
	private:
		double sx;
		double sy;
		double amp;
		double width;
};

/* Coulomb-like impurities in 3D space with a tunable exponent. */
class CoulombImpurities : public Noise {
	public:
		typedef std::tr1::tuple<double, double, double, double> impurity; // x-coord, y-coord, z-coord, alpha
		CoulombImpurities(double d, double e, double a, double maxd, DataLayout const& dl, Constraint const& constr, RNG& rng);
		void add_noise(DataLayout const& dl, double* pot_values) const;
		void write_realization_data(std::vector<double>& vec) const;
		DataLayout const& datalayout;
		Constraint const& constraint;
	private:
		double density; // density of impurities in three dimensions
		double exponent;
		double alpha;	// the potential of a single impurity is alpha/r^exp
		double maxd;	// the impurities are located furthest at maxd distance units from the calculation plane
		void init(RNG& rng);
		std::vector<impurity> impurities;
};

// TODO: There is a lot of code duplication going on in GaussianNoise,
// CoulombImpurities and HemisphereImpurities. Replace them with a common
// UniformImpurities with clever template code to set the type of a single
// impurity

/* Hemisphere-shaped finite-range impurities. */
class HemisphereImpurities : public Noise {
	public:
		typedef std::tr1::tuple<double, double, double, double> impurity; // x-coord, y-coord, amplitude, radius
		HemisphereImpurities(double d, double A, double r, DataLayout const& dl, Constraint const& constr, RNG& rng);
		void add_noise(DataLayout const& dl, double* pot_values) const;
		void write_realization_data(std::vector<double>& vec) const;
		DataLayout const& datalayout;
		Constraint const& constraint;
	private:
		double density; // density of impurities in two dimensions
		double amplitude;
		double radius;
		void init(RNG& rng);
		std::vector<impurity> impurities;
};

#endif // _NOISE_HPP_
