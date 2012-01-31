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

#include "datalayout.hpp"

DataLayout::DataLayout(size_t sx, size_t sy, double _dx) :
		sizex(sx), sizey(sy), N(sx*sy), dx(_dx),
		lenx(static_cast<double>(sx)*dx),
		leny(static_cast<double>(sy)*dx),
		posx(new double[sx]),
		posy(new double[sy]) {
	const int isx = static_cast<int>(sx);
	const int isy = static_cast<int>(sy);
	// The position coordinates will refer to the centers of each grid box
	for (int x=0; x < isx; x++) {
		posx[x] = static_cast<double>(2*x + 1 - isx)*0.5*dx;
	}
	for (int y=0; y < isy; y++) {
		posy[y] = static_cast<double>(2*y + 1 - isy)*0.5*dx;
	}
}

DataLayout::~DataLayout() {
	delete[] posx;
	delete[] posy;
}
