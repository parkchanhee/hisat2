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

bool isConcordant(long long int &location1, bool &forward1, long long int &location2, bool &forward2) {
    if (abs(location1-location2) > 500000) { return false; }
    if (location1 < location2) {
        if (forward1 && !forward2) { return true; }
    } else {
        if (!forward1 && forward2) { return true; }
    }
    return false;
}

int calculatePairScore_DNA (long long int &location0, int& AS0, bool& forward0, long long int &location1, int &AS1, bool &forward1, bool& concordant) {
    // this is the basic function to calculate pair score.
    // if the distance between 2 alignment is more than 1000, we reduce the score by the distance/100.
    // if two alignment is concordant we add 500000 to make sure to select the concordant pair as best pair.
    int score = 100*AS0 + 100*AS1;
    int distance = abs(location0 - location1);
    if (distance > 500000) { return numeric_limits<int>::min(); }
    if (distance > 1000) { score -= distance/100; }
    concordant = isConcordant(location0, forward0, location1, forward1);
    if (concordant) { score += 500000; }
    return score;
}

int calculatePairScore_RNA (long long int &location0, int& XM0, bool& forward0, long long int &location1, int &XM1, bool &forward1, bool& concordant) {
    // this is the basic function to calculate pair score.
    // if the distance between 2 alignment is more than 100,000, we reduce the score by the distance/1000.
    // if two alignment is concordant we add 500,000 to make sure to select the concordant pair as best pair.
    int score = -100*XM0 + -100*XM1;
    int distance = abs(location0 - location1);
    if (distance > 500000) { return numeric_limits<int>::min(); }
    if (distance > 100000) { score -= distance/1000; }
    concordant = isConcordant(location0, forward0, location1, forward1);
    if (concordant) { score += 500000; }
    return score;
}

bool MappingPosition::operator==(Alignment* o) {

    BTString* testChromosome;
    if (!o->repeat && o->pairToRepeat) {
        testChromosome = &o->pairToChromosome;
    } else {
        testChromosome = &o->chromosomeName;
    }
    //BTString* chromosome = (!o->repeat && o->pairToRepeat)?&o->pairToChromosome:&o->chromosomeName
    return (*locations[o->pairSegment] == o->location) &&
        (*locations[1-(o->pairSegment)] == o->pairToLocation) &&
        (*chromosome == *testChromosome);
    /*if (o->pairSegment == 0) {
        return (locations[0] == o->location) && (locations[1] == o->pairToLocation) && (chromosome == o->chromosomeName);
    } else {
        return (locations[0] == o->pairToLocation) && (locations[1] == o->location) && (chromosome == o->chromosomeName);
    }*/
}

/*MappingPosition::MappingPosition(Alignment* newAlignment) {
    locations[newAlignment->pairSegment] = &newAlignment->location;
    locations[1-newAlignment->pairSegment] = &newAlignment->pairToLocation;
    segmentExist[newAlignment->pairSegment] = true;
    chromosome = &newAlignment->chromosomeName;
    *//*if (newAlignment->repeat || newAlignment->pairToRepeat) {
        repeat = true;
    }*//*
    pairScore = numeric_limits<int>::min();
}

MappingPosition::MappingPosition (RepeatMappingPosition* repeat0, Alignment* newAlignment0, RepeatMappingPosition* repeat1=NULL, Alignment* newAlignment1=NULL) {
    locations[newAlignment0->pairSegment] = &repeat0->repeatLocation;
    chromosome = &repeat0->repeatChromosome;
    repeats[0] = repeat0;
    repeats[1] = repeat1;
    alignments[0] = newAlignment0;
    alignments[1] = newAlignment1;
    segmentExist[0] = true;
    if (alignments[1] != NULL) {
        locations[newAlignment1->pairSegment] = &repeat1->repeatLocation;
        segmentExist[1] = true;
    }
    AS = repeat0->AS;
    repeat = true;
}*/

void MappingPosition::install(Alignment *newAlignment) {
    locations[newAlignment->pairSegment] = &newAlignment->location;
    locations[1-newAlignment->pairSegment] = &newAlignment->pairToLocation;
    segmentExist[newAlignment->pairSegment] = true;
    chromosome = &newAlignment->chromosomeName;
    /*if (newAlignment->repeat || newAlignment->pairToRepeat) {
        repeat = true;
    }*/
    pairScore = numeric_limits<int>::min();
}

void MappingPosition::install (RepeatMappingPosition* repeat0, Alignment* newAlignment0, RepeatMappingPosition* repeat1=NULL, Alignment* newAlignment1=NULL) {
    locations[newAlignment0->pairSegment] = &repeat0->repeatLocation;
    chromosome = &repeat0->repeatChromosome;
    repeats[0] = repeat0;
    repeats[1] = repeat1;
    alignments[0] = newAlignment0;
    alignments[1] = newAlignment1;
    segmentExist[0] = true;
    if (alignments[1] != NULL) {
        locations[newAlignment1->pairSegment] = &repeat1->repeatLocation;
        segmentExist[1] = true;
    }
    AS = repeat0->AS;
    repeat = true;
}

bool MappingPositions::positionExist_new (Alignment* newAlignment) {
    if (positions.empty()) {
        index = 0;
        return false;
    }

    if (*positions[index] == newAlignment) {
        return positions[index]->segmentExist[newAlignment->pairSegment];
    }

    /*long long int* location0;
    long long int* location1;*/
    int segment = newAlignment->pairSegment;
    long long int* targetLocations[2];
    targetLocations[segment] = &newAlignment->location;
    targetLocations[1-segment] = &newAlignment->pairToLocation;

    /*if (newAlignment->pairSegment == 0) {
        location0 = &newAlignment->location;
        location1 = &newAlignment->pairToLocation;
    } else {
        location0 = &newAlignment->pairToLocation;
        location1 = &newAlignment->location;
    }*/

    return findPosition_new(targetLocations,
                            (!newAlignment->repeat && newAlignment->pairToRepeat)?newAlignment->pairToChromosome:newAlignment->chromosomeName,
                            segment);

}

/*bool MappingPositions::positionExist (Alignment* newAlignment, int& index) {
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
            !positions[i].segmentExist[*segment]) {
            index = i;
            return true;
        }
    }
    return false;
}*/

/*bool MappingPositions::positionExist (Alignment* newAlignment) {
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
            !positions[i].segmentExist[*segment]) {
            return true;
        }
    }
    return false;
}*/

/*bool MappingPositions::append (Alignment* newAlignment) {
    // return true if the position is not exist and will append to positions, else return false.

    int index;
    if (positionExist(newAlignment, index)) {
        return false;
    } else {
        //positions.push_back(MappingPosition(newAlignment->location, newAlignment->pairToLocation, newAlignment->chromosomeName, newAlignment->pairSegment));
        positions.emplace_back(newAlignment->location, newAlignment->pairToLocation, newAlignment->chromosomeName, newAlignment->pairSegment);
        return true;
    }
}*/

bool MappingPositions::append(Alignment* newAlignment) {

    if (positionExist_new(newAlignment)) {
        return false;
    } else {
        int segment = newAlignment->pairSegment;
        if (!positions.empty() && *positions[index] == newAlignment) {
            positions[index]->segmentExist[segment] = true;
            if (positions[index]->badAlignment) {
                return false;
            }
            positions[index]->alignments[segment] = newAlignment;
        } else {
            BTString* inputChromosome;
            if (!newAlignment->repeat && newAlignment->pairToRepeat) {
                inputChromosome = &newAlignment->pairToChromosome;
            } else {
                inputChromosome = &newAlignment->chromosomeName;
            }
            appendPosition(newAlignment->location, newAlignment->pairToLocation, inputChromosome, newAlignment->pairSegment);
            //positions.emplace_back(newAlignment->location, newAlignment->pairToLocation, inputChromosome, newAlignment->pairSegment);
            index = positions.size()-1;
            positions[index]->alignments[segment] = newAlignment;
            if (oppositeAlignment != NULL) {
                positions[index]->alignments[1-segment] = oppositeAlignment;
                positions[index]->segmentExist[1-segment] = true;
            }
        }
        return true;
    }
}

/*void MappingPositions::directAppend(Alignment *newAlignment, int& index) {
    if (index == -1) {
        positions.emplace_back(newAlignment->location, newAlignment->pairToLocation, newAlignment->chromosomeName, newAlignment->pairSegment);
    } else {
        positions[index].segmentExist[newAlignment->pairSegment] = true;
    }
}*/

void MappingPositions::outputPair(BTString& o) {
    int outputCount = 0;
    bool primary = true;
    for (int i = 0; i < positions.size(); i++) {
        if (positions[i]->pairScore == bestPairScore) {
            outputCount++;
            assert(positions[i]->alignments[0] != NULL);
            assert(positions[i]->alignments[1] != NULL);
            bool concordant = isConcordant(*positions[i]->locations[0],
                                           positions[i]->alignments[0]->forward,
                                           *positions[i]->locations[1],
                                           positions[i]->alignments[1]->forward);
            positions[i]->alignments[0]->setConcordant(concordant);
            positions[i]->alignments[1]->setConcordant(concordant);
            if (!positions[i]->repeat) {
                if (positions[i]->alignments[0]->outputted && positions[i]->alignments[1]->outputted) {
                    positions[i]->alignments[1]->outputted = false;
                }
                positions[i]->alignments[0]->setYS(positions[i]->alignments[1]);
                positions[i]->alignments[1]->setYS(positions[i]->alignments[0]);
                positions[i]->alignments[0]->outputRegularAlginemnt(o, positions[i]->locations[1], primary);
                positions[i]->alignments[1]->outputRegularAlginemnt(o, positions[i]->locations[0], primary);
            } else {
                //output repeat
                if (positions[i]->repeats[0]->outputted && positions[i]->repeats[1]->outputted) {
                    positions[i]->repeats[1]->outputted = false;
                }
                positions[i]->repeats[0]->setYS(positions[i]->repeats[1]);
                positions[i]->repeats[1]->setYS(positions[i]->repeats[0]);
                positions[i]->alignments[0]->outputRepeatAlignment(o, positions[i]->repeats[0], positions[i]->locations[1], primary);
                positions[i]->alignments[1]->outputRepeatAlignment(o, positions[i]->repeats[1], positions[i]->locations[0], primary);
            }
            primary = false;
        }
    }
    assert(outputCount == nBestPair);
}

void MappingPositions::outputSingle(BTString &o) {
    int outputCount = 0;
    bool primary = true;
    for (int i = 0; i < positions.size(); i++) {
        if (positions[i]->AS == bestAS && !positions[i]->badAlignment) {
            outputCount++;
            assert(positions[i]->alignments[0] != NULL);
            if (!positions[i]->repeat) {
                positions[i]->alignments[0]->outputRegularAlginemnt(o, NULL, primary);
            } else {
                positions[i]->alignments[0]->outputRepeatAlignment(o, positions[i]->repeats[0], NULL, primary);
            }
            primary = false;
        }
    }
}

bool MappingPositions::updatePairScore_regular() {
    int nPair;
    int score;
    score = positions[index]->alignments[0]->calculatePairScore(positions[index]->alignments[1], nPair);
    if (score > bestPairScore) {
        bestPairScore = score;
        nBestPair = nPair;
        concordantExist = positions[index]->alignments[0]->concordant;
    } else if (score == bestPairScore) {
        nBestPair += nPair;
    } else { // the newPair Score is less than bestPairScore, label it
        badAligned();
        return false;
    }
    positions[index]->pairScore = score;
    return true;
}

bool MappingPositions::updateAS_regular() {
    if (isBad()) { return false; }
    if (!positions[index]->alignments[0]->mapped) { return true; }
    int AS = positions[index]->alignments[0]->AS;
    if (AS > bestAS) {
        bestAS = AS;
        nBestSingle = 1;
    } else if (AS == bestAS) {
        nBestSingle++;
    } else {
        badAligned();
        return false;
    }
    positions[index]->AS = AS;
    return true;
}

bool MappingPositions::updatePairScore_repeat() {
    Alignment* alignments[2];
    alignments[0] = positions[index]->alignments[0];
    alignments[1] = positions[index]->alignments[1];
    if ((!alignments[0]->mapped || !alignments[1]->mapped) &&
        (bestPairScore >= (numeric_limits<int>::min()/2 - 1))) {
        badAligned();
        return false;
    }
    RepeatMappingPosition *repeatPosition0;
    RepeatMappingPosition *repeatPosition1;
    RepeatMappingPosition *repeatFlag0;
    RepeatMappingPosition *repeatFlag1;
    bool forward[2];
    forward[0] = alignments[0]->forward;
    forward[1] = alignments[1]->forward;
    bool DNA = alignments[0]->DNA;
    int score;
    bool concordant;
    for (int i = 0; i < alignments[0]->repeatPositions.size(); i++) {
        repeatPosition0 = &alignments[0]->repeatPositions.positions[i];
        repeatFlag0 = repeatPosition0->repeatFlagInfo==NULL ? repeatPosition0 : repeatPosition0->repeatFlagInfo;
        for (int j = 0; j < alignments[1]->repeatPositions.size(); j++) {
            repeatPosition1 = &alignments[1]->repeatPositions.positions[j];
            if (repeatPosition0->repeatChromosome == repeatPosition1->repeatChromosome) {
                repeatFlag1 = repeatPosition1->repeatFlagInfo==NULL ? repeatPosition1 : repeatPosition1->repeatFlagInfo;
                if (DNA) {
                    score = calculatePairScore_DNA(repeatPosition0->repeatLocation,
                                                   repeatFlag0->AS,
                                                   forward[0],
                                                   repeatPosition1->repeatLocation,
                                                   repeatFlag1->AS,
                                                   forward[1],
                                                   concordant);
                } else {
                    score = calculatePairScore_RNA(repeatPosition0->repeatLocation,
                                                   repeatFlag0->XM,
                                                   forward[0],
                                                   repeatPosition1->repeatLocation,
                                                   repeatFlag1->XM,
                                                   forward[1],
                                                   concordant);
                }
                if (score >= bestPairScore) {
                    //positions.emplace_back(repeatPosition0, alignments[0], repeatPosition1, alignments[1]);
                    appendPosition(repeatPosition0, alignments[0], repeatPosition1, alignments[1]);
                    positions.back()->pairScore = score;
                    if (score > bestPairScore) {
                        nBestPair = 1;
                        bestPairScore = score;
                        concordantExist = concordant;
                    } else {
                        nBestPair++;
                    }
                }
            }
        }
    }
    return true;
}

bool MappingPositions::updateAS_repeat() {
    if (isBad()) { return false; }
    Alignment* alignment = positions[index]->alignments[0];
    RepeatMappingPosition* repeatPosition;
    badAligned(); // label this as bad alignment to avoid directly output.
    int AS;
    for (int i = 0; i < alignment->repeatPositions.size(); i++) {
        repeatPosition = &alignment->repeatPositions.positions[i];
        AS = (repeatPosition->repeatFlagInfo == NULL)?repeatPosition->AS : repeatPosition->repeatFlagInfo->AS;
        if (AS >= bestAS) {
            appendPosition(repeatPosition, alignment);
            //positions.emplace_back(repeatPosition, alignment);
            //positions.back().AS = AS;
            if (AS > bestAS) {
                bestAS = AS;
                nBestSingle = 1;
            } else {
                nBestSingle++;
            }
        }
    }
    return true;
}


bool MappingPositions::updatePairScore() {
    if (!mateExist()) { return true; }

    assert(positions[index]->alignments[0] != NULL);
    assert(positions[index]->alignments[1] != NULL);

    if (positions[index]->alignments[0]->repeat || positions[index]->alignments[1]->repeat) {
        return updatePairScore_repeat();
    } else {
        return updatePairScore_regular();
    }
}

bool MappingPositions::updateAS() {
    if (positions[index]->alignments[0]->repeat) {
        return updateAS_repeat();
    } else {
        return updateAS_regular();
    }
}

void MappingPositions::appendPosition(Alignment *newAlignment) {
    MappingPosition* newPostion;
    getFreePositionPointer(newPostion);
    newPostion->install(newAlignment);
    positions.push_back(newPostion);
}

void MappingPositions::appendPosition(RepeatMappingPosition* repeat0, Alignment* newAlignment0, RepeatMappingPosition* repeat1, Alignment* newAlignment1) {
    MappingPosition* newPostion;
    getFreePositionPointer(newPostion);
    newPostion->install(repeat0, newAlignment0, repeat1, newAlignment1);
    positions.push_back(newPostion);
}

int Alignment::calculatePairScore(Alignment *inputAlignment, int &nPair) {
    // calculate the pairScore for a pair of alignment result. Output pair Score and number of pair (1 for non-repeat).
    // Do not update their pairScore.

    int pairScore = numeric_limits<int>::min();
    nPair = 0;
    bool concordant;
    if (pairSegment == inputAlignment->pairSegment){
        // when 2 alignment results are from same pair segment, output the lowest score and number of pair equal zero.
        pairScore = numeric_limits<int>::min();
    } else if (!mapped && !inputAlignment->mapped) {
        // both ummaped.
        pairScore = numeric_limits<int>::min()/2 - 1;
    } else if (!mapped || !inputAlignment->mapped) {
        // one of the segment unmapped.
        pairScore = numeric_limits<int>::min()/2;
        nPair = 1;
    } else if ((!repeat && !inputAlignment->repeat)){
        // both mapped and (both non-repeat or not expand repeat)
        bool concordant;
        if (DNA) {
            pairScore = calculatePairScore_DNA(location,
                                               AS,
                                               forward,
                                               inputAlignment->location,
                                               inputAlignment->AS,
                                               inputAlignment->forward,
                                               concordant);
        } else {
            pairScore = calculatePairScore_RNA(location,
                                               XM,
                                               forward,
                                               inputAlignment->location,
                                               inputAlignment->XM,
                                               inputAlignment->forward,
                                               concordant);
        }

        setConcordant(concordant);
        inputAlignment->setConcordant(concordant);

        nPair = 1;
    }
    return pairScore;
}


