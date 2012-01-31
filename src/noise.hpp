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

#include <vector>
#include <tr1/tuple>

#include "itp2d_common.hpp"
#include "datalayout.hpp"
#include "rng.hpp"

/* Interface class representing "noise" which can be added to some array of
 * values, such as a potential. */
class Noise {
	public:
		virtual ~Noise() {};
		virtual void add_noise(DataLayout const& dl, double* pot_values, RNG& rng) const = 0;
		virtual std::string const& get_description() const { return description; }
	protected:
		std::string description;
};

/* The simple case of no noise at all */
class NoNoise : public Noise {
	public:
		NoNoise() { description = "none"; }
		void add_noise(__attribute__((unused)) DataLayout const& dl, __attribute__((unused))double* pot_values, __attribute__((unused))RNG& rng) const {}
};

/* Randomly distributed gaussian spikes with normally distributed width
 * and normally distributed amplitude. */
class GaussianSpikes : public Noise {
	public:
		GaussianSpikes(double density, double amplitude_mean, double amplitude_stdev, double width_mean, double width_stdev);
		void add_noise(DataLayout const& dl, double* pot_values, RNG& rng) const;
	private:
		const double density;
		const double amplitude_mean;
		const double amplitude_stdev;
		const double width_mean;
		const double width_stdev;
};

#endif // _NOISE_HPP_
