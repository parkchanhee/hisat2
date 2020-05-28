//
// Created by s429377 on 1/10/20.
//

#include "TLA.h"

bool MappingPositions::append (Alignment* newAlignment) {
    // return true if the position is not exist and will append to positions, else, false.
    // if alignment is repeat (mapped to repeat index), don't push to positions, return true.

    long long int location = newAlignment->location;
    string chromosome = newAlignment->chromosomeName.toZBuf();
    int pairSegment = newAlignment->pairSegment;
    int AS = newAlignment->AS;
    bool concordant = newAlignment->concordant;

    if (newAlignment->repeat) {
        return true;
    }

    int index;
    if (positionExist(location, chromosome, pairSegment, index)) {
        if ((positions[index].concordant == false) && (concordant == true)) {
            return true;
        } else {
            return false;
        }
    } else {
        positions.push_back(MappingPosition(location, chromosome, AS, pairSegment, concordant));
        return true;
    }
}

