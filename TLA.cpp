/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "TLA.h"

bool MappingPositions::append (Alignment* newAlignment) {
    // return true if the position is not exist and will append to positions, else, false.
    // if alignment is repeat (mapped to repeat index), don't push to positions, return true.

    long long int location = newAlignment->location;
    string chromosome = newAlignment->chromosomeName.toZBuf();
    int pairSegment = newAlignment->pairSegment;
    bool concordant = newAlignment->concordant;

    /*if (newAlignment->repeat) {
        return true;
    }*/

    int index;
    if (positionExist(location, chromosome, pairSegment, index)) {
        return (!positions[index].concordant) && concordant;
    } else {
        positions.push_back(MappingPosition(location, chromosome, pairSegment, concordant));
        return true;
    }
}

bool positionExist (Alignment* newAlignment) {

}

