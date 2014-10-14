/* Copyright 2014 Perttu Luukko

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
 * Unit tests for the Timer class.
 */

#include "test_timer.hpp"

// Test that the timer produces sensible results
TEST(timer, basic_timing) {
	const double tolerance = 1e-3;
	const unsigned int timings = 5;
	const unsigned int sleep_us = 10000;
	const double expected_time = timings*sleep_us*1e-6;
	Timer timer;
	for (unsigned i=0; i<timings; i++) {
		timer.start();
		usleep(sleep_us);
		timer.stop();
	}
	const double elapsed_time = timer.get_time();
	EXPECT_NEAR(expected_time, elapsed_time, tolerance);
}
