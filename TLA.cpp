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

bool MappingPositions::positionExist (Alignment* newAlignment, int& index) {
    for (int i = 0; i < positions.size(); i++) {
        if ((positions[i].location == newAlignment->location) &&
            (positions[i].chromosome == newAlignment->chromosomeName) &&
            (positions[i].pairSegment == newAlignment->pairSegment)) {
            index = i;
            return true;
        }
    }
    return false;
}

bool MappingPositions::positionExist (Alignment* newAlignment) {
    for (int i = 0; i < positions.size(); i++) {
        if ((positions[i].location == newAlignment->location) &&
            (positions[i].chromosome == newAlignment->chromosomeName) &&
            (positions[i].pairSegment == newAlignment->pairSegment)) {
            return (!positions[i].concordant) && newAlignment->concordant;
        }
    }
    return false;
}

bool MappingPositions::append (Alignment* newAlignment) {
    // return true if the position is not exist and will append to positions, else return false.

    int index;
    if (positionExist(newAlignment, index)) {
        return (!positions[index].concordant) && newAlignment->concordant;
    } else {
        positions.push_back(MappingPosition(newAlignment->location, newAlignment->chromosomeName, newAlignment->pairSegment, newAlignment->concordant));
        return true;
    }
}

void MappingPositions::directAppend(Alignment *newAlignment) {
    positions.push_back(MappingPosition(newAlignment->location, newAlignment->chromosomeName, newAlignment->pairSegment, newAlignment->concordant));
}



