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
    if (positions.empty()) {
        return false;
    }

    long long int* location0;
    long long int* location1;
    int* segment = &newAlignment->pairSegment;

    if (newAlignment->pairSegment == 0) {
        location0 = &newAlignment->location;
        location1 = &newAlignment->pairToLocation;
    } else {
        location0 = &newAlignment->pairToLocation;
        location1 = &newAlignment->location;
    }

    for (int i = 0; i < positions.size(); i++) {
        if ((positions[i].location0 == *location0) &&
            (positions[i].location1 == *location1) &&
            (positions[i].chromosome == newAlignment->chromosomeName) &&
            positions[i].segmentExist[*segment]) {
            index = i;
            return true;
        }
    }
    return false;
}

bool MappingPositions::positionExist (Alignment* newAlignment) {
    if (positions.empty()) {
        return false;
    }

    long long int* location0;
    long long int* location1;
    int* segment = &newAlignment->pairSegment;

    if (*segment == 0) {
        location0 = &newAlignment->location;
        location1 = &newAlignment->pairToLocation;
    } else {
        location0 = &newAlignment->pairToLocation;
        location1 = &newAlignment->location;
    }

    for (int i = 0; i < positions.size(); i++) {
        if ((positions[i].location0 == *location0) &&
            (positions[i].location1 == *location1) &&
            (positions[i].chromosome == newAlignment->chromosomeName) &&
            positions[i].segmentExist[*segment]) {
            return true;
        }
    }
    return false;
}

bool MappingPositions::append (Alignment* newAlignment) {
    // return true if the position is not exist and will append to positions, else return false.

    int index;
    if (positionExist(newAlignment, index)) {
        return false;
    } else {
        //positions.push_back(MappingPosition(newAlignment->location, newAlignment->pairToLocation, newAlignment->chromosomeName, newAlignment->pairSegment));
        positions.emplace_back(newAlignment->location, newAlignment->pairToLocation, newAlignment->chromosomeName, newAlignment->pairSegment);
        return true;
    }
}

void MappingPositions::directAppend(Alignment *newAlignment, int& index) {
    if (index == -1) {
        positions.emplace_back(newAlignment->location, newAlignment->pairToLocation, newAlignment->chromosomeName, newAlignment->pairSegment);
    } else {
        positions[index].segmentExist[newAlignment->pairSegment] = true;
    }

}



