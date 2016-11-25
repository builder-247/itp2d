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
		virtual void new_realization(std::list<double>& params) const = 0;
		virtual void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const = 0;
		virtual size_t get_num_params() const = 0;
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
		virtual Constraint const& get_constraint() const = 0;
	protected:
		std::list<coordinate_pair> coordinates;
		std::string description;
};

// Regular spatial noise class built from impurity type + impurity distribution

class SpatialImpurities : public Noise {
	public:
		typedef ImpurityDistribution::coordinate_pair coordinate_pair;
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
			std::stringstream ss;
			ss << "Gaussian spikes, normally distributed amplitude with mean " << amp_mean
				<< " and standard deviation " << amp_stdev
				<< ", normally distributed width with mean " << width_mean
				<< " and standard deviation " << width_stdev;
			description = ss.str();
		}
		void new_realization(std::list<double>& params) const {
			const double A = amp_mean + amp_stdev*rng.gaussian_rand();
			const double w = width_mean + width_stdev*rng.gaussian_rand();
			params.push_back(A);
			params.push_back(w);
		}
		void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const;
		size_t get_num_params() const { return num_params; }
	private:
		double amp_mean;
		double amp_stdev;
		double width_mean;
		double width_stdev;
		RNG& rng;
};

class CoulombImpurities : public ImpurityType {
	static const size_t num_params = 0;
	public:
		CoulombImpurities(double e, double _alpha, RNG& _rng) :
				halfexponent(e/2), alpha(_alpha), rng(_rng) {
			std::stringstream ss;
			ss << "Coulomb-like impurities with exponent " << 2*halfexponent
				<< " and prefactor " << alpha;
			description = ss.str();
		}
		void new_realization(__attribute__((unused)) std::list<double>& params) const {}
		void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const;
		size_t get_num_params() const { return num_params; }
	private:
		double halfexponent;
		double alpha;
		RNG& rng;
};

class HemisphereImpurities : public ImpurityType {
	static const size_t num_params = 0;
	public:
		HemisphereImpurities(double amplitude, double radius, RNG& _rng) :
				A(amplitude), r2(radius*radius), rng(_rng) {
			std::stringstream ss;
			ss << "hemisphere bumps with amplitude " << A
				<< " and radius " << std::sqrt(r2);
			description = ss.str();
		}
		void new_realization(__attribute__((unused)) std::list<double>& params) const {}
		void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const;
		size_t get_num_params() const { return num_params; }
	private:
		double A;
		double r2;
		RNG& rng;
};

class DeltaImpurities : public ImpurityType {
	static const size_t num_params = 1;
	public:
		DeltaImpurities(double _amp_mean, double _amp_stdev, DataLayout const& dl, RNG& _rng) :
				amp_mean(_amp_mean), amp_stdev(_amp_stdev), inverse_unit_area(1.0/(dl.dx*dl.dx)), rng(_rng) {
			std::stringstream ss;
			ss << "Delta spikes, normally distributed integral with mean " << amp_mean
				<< " and standard deviation " << amp_stdev;
			description = ss.str();
		}
		void new_realization(std::list<double>& params) const {
			const double A = amp_mean + amp_stdev*rng.gaussian_rand();
			params.push_back(A);
		}
		void add_noise(double x, double y, std::vector<double> const& params, DataLayout const& dl, double* pot_values) const;
		size_t get_num_params() const { return num_params; }
	private:
		double amp_mean;
		double amp_stdev;
		double inverse_unit_area;
		RNG& rng;
};

// Individual distribution types

class UniformImpurities : public ImpurityDistribution {
	public:
		UniformImpurities(double density, DataLayout const& dl, Constraint const& _constraint, RNG& rng) : 
				constraint(_constraint) {
			std::stringstream ss;
			ss << "uniform distribution with density " << density << " impurities per unit square";
			description = ss.str();
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


class FixedImpurities : public ImpurityDistribution {
	public:
		FixedImpurities(std::list<coordinate_pair> _coordinates) {
			description = "user-defined locations of impurities";
			coordinates = _coordinates;
		}
		Constraint const& get_constraint() const { return constraint; }
	private:
		NoConstraint constraint;
};

// Parser functions for returning instances from user-provided desciption
// strings

ImpurityType const* parse_impurity_type_description(std::string const& type_str, DataLayout const& dl, RNG& rng);

ImpurityDistribution const* parse_impurity_distribution_description(std::string const& distribution_str, DataLayout const& dl, Constraint const& constraint, RNG& rng);

#endif // _NOISE_HPP_
