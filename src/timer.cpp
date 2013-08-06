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

#include "timer.hpp"

Timer::Timer() : elapsed_nsec(0), running(false) {
	// Test that CLOCK_MONOTONIC works
	timespec* t = new timespec;
	const int retval = clock_gettime(CLOCK_MONOTONIC, t);
	if (retval != 0)
		throw GeneralError("CLOCK_MONOTONIC not supported. Cannot get reliable timing data");
	delete t;
}
