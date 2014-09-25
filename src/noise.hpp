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
		virtual std::string const& get_distribution_description() const = 0;
		virtual std::string const& get_constraint_description() const = 0;
		// Write internal data which can be used to re-create the noise realization. For example for
		// Gaussian impurities, store the positions, amplitudes and widths of the Gaussians.
		virtual void write_realization_data(std::vector<double>& vec) const = 0;
};

// Interface class of impurity types

class ImpurityType {
	public:
		virtual ~ImpurityType() {};
		virtual void new_realization(std::vector<double>& params) = 0;
		virtual void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const = 0;
		std::string const& get_description() const { return description; }
	protected:
		std::string description;
};

// Interface class of impurity distributions

class ImpurityDistribution {
	public:
		typedef std::pair<double, double> coordinate_pair;
		virtual ~ImpurityDistribution() {};
		std::list<coordinate_pair> const& get_coordinates() const { return coordinates; }
		std::string const& get_description() const { return description; }
		virtual Constraint const& get_constraint() const { return NoConstraint(); }
	protected:
		std::list<coordinate_pair> coordinates;
		std::string description;
};

// Regular spatial noise class built from impurity type + impurity distribution

class SpatialImpurities : public Noise {
	public:
		SpatialImpurities(ImpurityType const& type, ImpurityDistribution const& distribution);
		void add_noise(DataLayout const& dl, double* pot_values) const;
		std::string const& get_description() const { return type.get_description(); }
		std::string const& get_distribution_description() const { return distribution.get_description(); }
		std::string const& get_constraint_description() const { return distribution.get_constraint().get_description(); }
		void write_realization_data(std::vector<double>& vec) const;
		ImpurityType const& type;
		ImpurityDistribution const& distribution;
	private:
		std::list<double> realization_data;
};

// The trivial case of no noise at all

class NoNoise : public Noise {
	public:
		NoNoise() {
			description = "none";
		}
		std::string const& get_description() const { return description; }
		std::string const& get_distribution_description() const { return description; }
		std::string const& get_constraint_description() const { return description; }
		void add_noise(__attribute__((unused)) DataLayout const& dl, __attribute__((unused)) double* pot_values) const {}
		void write_realization_data(std::vector<double>& vec) const { vec.clear(); }
	private:
		std::string description;
};

// Individual impurity types

class GaussianImpurities : public ImpurityType {
	static const size_t num_params = 2;
	public:
		GaussianImpurities(double _amp_mean, double _amp_stdev, double _width_mean, double _width_stdev, RNG& _rng) :
				amp_mean(_amp_mean), amp_stdev(_amp_stdev), width_mean(_width_mean), width_stdev(_width_stdev), rng(_rng) {
			description = "TODO";
		}
		void new_realization(std::vector<double>& params) {
			const double A = amp_mean + amp_stdev*rng.gaussian_rand();
			const double w = width_mean + width_stdev*rng.gaussian_rand();
			params.push_back(A);
			params.push_back(w);
		}
		void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const;
	private:
		double amp_mean;
		double amp_stdev;
		double width_mean;
		double width_stdev;
		RNG& rng;
};

// Individual distribution types

class UniformImpurities : public ImpurityDistribution {
	public:
		UniformImpurities(double density, DataLayout const& dl, Constraint const& _constraint, RNG& rng) : 
				constraint(_constraint) {
			description = "TODO";
			// Draw the number of impurities from a Poisson distribution
			const double lambda = density*dl.lenx*dl.leny;
			const unsigned int N = rng.poisson_rand(lambda);
			// Randomly distribute the positions of the impurities
			for (unsigned int n=0; n<N; n++) {
				const double x = (rng.uniform_rand()-0.5)*dl.lenx;
				const double y = (rng.uniform_rand()-0.5)*dl.leny;
				if (constraint.check(x, y)) {
					coordinates.push_back(std::make_pair(x, y));
				}
			}
		}
		Constraint const& get_constraint() const { return constraint; }
		Constraint const& constraint;
};


// A parser functions for returning instances from a user provided desciption
// strings

ImpurityType const* parse_impurity_type_description(std::string const& type_str, DataLayout const& dl, RNG& rng);

ImpurityDistribution const* parse_impurity_distribution_description(std::string const& distribution_str, DataLayout const& dl, Constraint const& constraint, RNG& rng);

// Old noise types
/*
* Randomly distributed gaussian spikes with normally distributed width
 * and normally distributed amplitude.
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

// Coulomb-like impurities in 3D space with a tunable exponent.
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

// Hemisphere-shaped finite-range impurities.
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
*/

#endif // _NOISE_HPP_
