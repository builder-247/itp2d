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
 * A simple wrapper class for random number generation based on TR1/random.
 */

#ifndef _RNG_HPP_
#define _RNG_HPP_

#ifndef _BSD_SOURCE
#define _BSD_SOURCE
#endif
#include <tr1/random>
#include <sys/time.h>

class RNG {
public:
	RNG();	// If no seed is given, one is generated based on system time.
	RNG(unsigned long int seed);
	inline double gaussian_rand() { return normal_distribution(uniform_rng); }
	inline double uniform_rand() { return uniform_rng(); }
	inline bool bernoulli_trial(double p) { return bernoulli_distribution_type(p)(uniform_rng); }
	inline unsigned long int get_seed() const { return seed; }
	typedef std::tr1::mt19937 base_rng_type;
	typedef std::tr1::normal_distribution<double> normal_distribution_type;
	typedef std::tr1::uniform_real<double> uniform_distribution_type;
	typedef std::tr1::bernoulli_distribution bernoulli_distribution_type;
	typedef std::tr1::variate_generator<base_rng_type, uniform_distribution_type> uniform_rng_type;
	static unsigned long int produce_random_seed();
private:
	unsigned long int seed;
	// We must have a uniform double-valued generator not only to generate
	// uniform random numbers, but also to generate normal distributed random
	// numbers. This is because due to limitations in GCC's implementation of
	// TR1 we cannot just call normal_distribution(base_rng), because the base
	// mt19937 is integer valued and those don't work with some distributions.
	// So we first create a double-valued generator with variate_generator, and
	// use that as a backend to get normal distributed values.
	uniform_rng_type uniform_rng;
	normal_distribution_type normal_distribution;
};

#endif // _RNG_HPP_
