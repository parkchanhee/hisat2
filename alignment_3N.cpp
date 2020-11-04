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

#include "alignment_3N.h"
#include "aln_sink.h"

bool Alignment::isConcordant(long long int &location1, bool &forward1, long long int &location2, bool &forward2) {
    if (abs(location1-location2) > 500000) { return false; }
    if (location1 < location2) {
        if (forward1 && !forward2) { return true; }
    } else {
        if (!forward1 && forward2) { return true; }
    }
    return false;
}

int Alignment::calculatePairScore_DNA (long long int &location0, int& AS0, bool& forward0, long long int &location1, int &AS1, bool &forward1, bool& concordant) {
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

int Alignment::calculatePairScore_RNA (long long int &location0, int& XM0, bool& forward0, long long int &location1, int &XM1, bool &forward1, bool& concordant) {
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

/**
 * calculate the pairScore for a pair of alignment result. Output pair Score and number of pair.
 * Do not update their pairScore.
 */
int Alignment::calculatePairScore(Alignment *inputAlignment, int &nPair) {
    int pairScore = numeric_limits<int>::min();
    nPair = 0;
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

void Alignments::reportStats_single(ReportingMetrics& met) {

    int nAlignment = alignmentPositions.nBestSingle;
    if (nAlignment == 0) {
        met.nunp_0++;
    } else {
        met.nunp_uni++;
        if (nAlignment == 1) { met.nunp_uni1++; }
        else { met.nunp_uni2++; }
    }
}

void Alignments::reportStats_paired(ReportingMetrics& met) {
    if (!alignmentPositions.concordantExist) {
        met.nconcord_0++;
        if (alignmentPositions.nBestPair == 0) {
            met.nunp_0_0 += 2;
            return;
        }
        if (alignmentPositions.bestPairScore == numeric_limits<int>::min()/2) {
            // one mate is unmapped, one mate is mapped
            met.nunp_0_0++;
            met.nunp_0_uni++;
            if (alignmentPositions.nBestPair == 1) { met.nunp_0_uni1++; }
            else { met.nunp_0_uni2++; }
        } else { //both mate is mapped
            if (alignmentPositions.nBestPair == 1) {
                met.ndiscord++;
                return;
            }
            else {
                met.nunp_0_uni += 2;
                met.nunp_0_uni2 += 2;
            }
        }
    } else {
        assert(alignmentPositions.nBestPair > 0);
        met.nconcord_uni++;
        if (alignmentPositions.nBestPair == 1) { met.nconcord_uni1++; }
        else { met.nconcord_uni2++; }
    }
}


