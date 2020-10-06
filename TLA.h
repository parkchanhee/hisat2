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

#ifndef HISAT2_TLA_H
#define HISAT2_TLA_H

#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <algorithm>
#include "sstring.h"
#include "util.h"
#include "hisat2lib/ht2.h"
#include "read.h"
#include "outq.h"
#include "reference.h"
#include <unistd.h>
#include <queue>

extern char convertedFrom;
extern char convertedTo;
extern char convertedFromComplement;
extern char convertedToComplement;
extern EList<ht2_handle_t> repeatHandles;
extern struct ht2_index_getrefnames_result *refNameMap;
extern int repeatLimit;
extern bool uniqueOutputOnly;

using namespace std;

class Alignment;

class ConvertMatrixTLA {
    char convertFrom = 'A';
    char convertTo = 'A';
    string allBase = "ACGT";
    string allBaseLower = "acgt";

    int charToInt(char inputChar) {
        return allBase.find(inputChar);
    }

    int complement(char inputChar) {
        return allBase[3-charToInt(inputChar)];
    }

    void convertMatrix() {
        restoreNormal();
        for (int i = 0; i < 4; i++) {
            char base = allBase[i];
            char lowerBase = allBaseLower[i];
            if (convertFrom == base) {
                asc2dna[base] = charToInt(convertTo);
                asc2dna[lowerBase] = charToInt(convertTo);
            } else if (complement(convertFrom) == base) {
                asc2dnacomp[base] = convertTo;
                asc2dnacomp[lowerBase] = convertTo;
                dnacomp[i] = charToInt(convertTo);
            }
        }
    }
public:
    ConvertMatrixTLA(){

    };

    void convert(char from, char to)  {
        convertFrom = from;
        convertTo = to;
        convertMatrix();
    }

    void inverseConversion() {
        convertFrom = complement(convertFrom);
        convertTo = complement(convertTo);
        convertMatrix();
    }

    void restoreNormal() {
        for (int i = 0; i < 4; i++) {
            char base = allBase[i];
            char lowerBase = allBaseLower[i];
            asc2dna[base] = charToInt(base);
            asc2dna[lowerBase] = charToInt(base);
            asc2dnacomp[base] = complement(base);
            asc2dnacomp[lowerBase] = complement(base);
            dnacomp[i] = charToInt(complement(base));
        }
    }

    void restoreConversion() {
        convertMatrix();
    }
};

class Cigar {
    int len;
    char label;
public:

    Cigar() {

    }

    Cigar(int inputLen, char inputLabel): len(inputLen), label(inputLabel) {
    }

    int& getLen() {
        return len;
    };

    char& getLabel() {
        return label;
    }
};

class RepeatPosition{
public:
    long long int location;
    BTString chromosome;
    int pairScore;
    long long int pairToLocation;
    bool concordant;
    BTString YZ;

    void initialize() {
        pairScore = numeric_limits<int>::min();
        pairToLocation = 0;
        concordant = false;
    }

    RepeatPosition() {
        initialize();
    }
    RepeatPosition(long long int &inputLocation, BTString &inputChromosome, BTString &repeatYZ):
        location(inputLocation), chromosome(inputChromosome),YZ(repeatYZ) {
        initialize();
    }
};

class MappingPosition {
public:
    long long int location;
    BTString chromosome;
    int AS;
    int pairScore;
    int pairSegment;
    long long int pairToLocation;
    bool concordant;

    // for repeat alignment only
    BTString MD;
    int XM;
    int NM;
    int YS;
    int TC;
    BTString MP;
    int RA_Array[5][5] = {{0,},};
    BTString refSequence;
    EList<RepeatPosition> repeatPositions;

    MappingPosition() {

    }

    MappingPosition (long long int inputLocation, BTString inputChromosome, int inputReadSegment=0, bool inputConcordant=false){
        concordant = inputConcordant;
        location = inputLocation;
        chromosome = inputChromosome;
        pairSegment = inputReadSegment;
        pairScore = numeric_limits<int>::min();
    }

    MappingPosition (long long int inputLocation, BTString inputChromosome, BTString &inputRefSequence, int &inputAS, BTString &inputMD, int &inputXM, int &inputNM, int &inputTC, BTString &inputMP, int inputRA[5][5], BTString &repeatYZ, int inputYS):
            location(inputLocation), chromosome(inputChromosome), refSequence(inputRefSequence), AS(inputAS), MD(inputMD), XM(inputXM), NM(inputNM), TC(inputTC), MP(inputMP), YS(inputYS){
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                RA_Array[i][j] = inputRA[i][j];
            }
        }
        pairScore = numeric_limits<int>::min();
        addRepeatPosition(inputLocation, inputChromosome, repeatYZ);
    }

    void addRepeatPosition(long long int &inputLocation, BTString &inputChromosome, BTString repeatYZ) {
        repeatPositions.push_back(RepeatPosition(inputLocation, inputChromosome, repeatYZ));
    }
};

class MappingPositions {
public:
    EList<MappingPosition> positions;
    int bestPairScore;
    int nPair;

    void initialize() {
        positions.clear();
        bestPairScore = numeric_limits<int>::min();
        nPair = 0;
    }

    MappingPositions() {
        initialize();
    };

    bool sequenceExist (BTString & refSequence, int &index) {
        // return true if reference sequence is exist, else, return false.
        for (int i = 0; i < positions.size(); i++) {
            if (refSequence == positions[i].refSequence) {
                index = i;
                return true;
            }
        }
        return false;
    }

    bool positionExist (Alignment* newAlignment, int& index);

    /*bool positionExist (long long int &location, string &chromosome, int pairSegment, int &index) {
        // return true if position is exist, else, return false.
        for (int i = 0; i < positions.size(); i++) {
            if ((positions[i].location == location) &&
                (positions[i].chromosome == chromosome) &&
                (positions[i].pairSegment == pairSegment)) {
                index = i;
                return true;
            }
        }
        return false;
    }*/

    void appendRepeat(BTString &chromosome, long long int &location, int &index) {
        positions[index].addRepeatPosition(location, chromosome, positions[index].repeatPositions[0].YZ);
    }

    void appendRepeat(long long int &location, BTString &chromosome, BTString &refSequence, int &AS, BTString &MD, int &XM, int &NM, int &TC, BTString &MP, int RA[5][5], BTString &repeatYZ, int YS=0) {
        positions.push_back(MappingPosition(location, chromosome, refSequence, AS, MD, XM, NM, TC, MP, RA, repeatYZ, YS));
    }

    bool append(Alignment* newAlignment);

    /*bool append (string chromosome, long long int location, int pairSegment=0) {
        // return true if the position is not exist and will append to positions, else, false.
        int index;
        if (positionExist(location, chromosome, pairSegment, index)) {
            return false;
        } else {
            positions.emplace_back(location, chromosome);
            return true;
        }
    }*/
};

class Alignment {
public:
    BTString readName;
    int flag;
    BTString chromosomeName;
    int chromosomeIndex;
    long long int location;
    BTString MAPQ;
    BTString cigarString;
    int cigarLength; // this is the length that read could cover genome.
    BTString pairToChromosome;
    long long int pairToLocation;
    long long int pairingDistance;
    BTString readSequence;
    BTString readQuality;
    EList<Cigar> cigarSegments;

    // save original sequence and quality score for output
    BTString originalFw;
    BTString readQualityFw;

    //tags
    int AS; // alignment score
    int NH; // number of alignment
    int XM; // number of mismatch
    int NM; // edit distance
    int YS; // mate's AS
    BTString MD;

    // extraTags
    int TC;
    BTString MP;
    int RA_Array[5][5] = {{0,},};
    BTString YZ;  // this tag shows alignment strand:
                // Original top strand (OT)
                // Complementary to original top strand (CTOT)
                // Complementary to original bottom strand (CTOB)
                // Original bottom strand (OB)

    // unChanged tags
    BTString unChangedTags;

    // intermediate variable
    //char* locationPointer; // point to the first char for mapping location in genome.
    //bool planA;
    int TLAcycle;
    bool paired;
    bool forward;
    bool primaryAlignment; // this is for output flag adjustment.
    bool mapped;
    bool concordant;
    int pairSegment; // 0 for first segment, 1 for second segment.
    struct ht2_repeat_expand_result *repeatResult = nullptr;
    int pairScore; // to identify the better pair

    EList<Alignment*> oppositePairAddresses;
    // for repeatAlignment only
    bool repeat;
    bool pairToRepeat;
    MappingPositions repeatPositions; // only have chromosme and location informations
    int bestAS;
    int nBestRepeat;

    int conversionCount[2] = {0};
    string intToBase = "ACGTN";

    void initialize() {
        readName.clear();
        flag = -1;
        chromosomeName.clear();
        location = 0;
        MAPQ.clear();
        cigarString.clear();
        cigarLength = 0;
        pairToChromosome.clear();
        pairToLocation = 0;
        pairingDistance = 0;
        readSequence.clear();
        readQuality.clear();
        chromosomeIndex = -1;
        cigarSegments.clear();

        AS = numeric_limits<int>::min();
        NH = 0;
        XM = 0;
        NM = 0;
        YS = 0;
        MD.clear();

        clearExtraTags();
        unChangedTags.clear();

        bestAS = numeric_limits<int>::min();
        repeat = false;
        pairToRepeat = false;
        nBestRepeat = 0;
        paired = false;
        pairScore = numeric_limits<int>::min();
        oppositePairAddresses.clear();
        repeatPositions.initialize();
        conversionCount[0] = 0;
        conversionCount[1] = 0;
        YZ.clear();
        TLAcycle = -1;

        if (repeatResult != nullptr) {
            free(repeatResult);
            repeatResult = nullptr;
        }
    }

    Alignment() {
        initialize();
    }

    ~Alignment() {
        if (repeatResult != nullptr) free(repeatResult);
    }

    void setConcordant(bool concordant_) {
        // change concordant status and flag
        concordant = concordant_;
        if ((flag&2) && !concordant) flag -= 2;
        else if (!(flag&2) && concordant) flag += 2;
    }
    void updateNH(int nAlignment) {
        if (!mapped) return;

        // update NH and MAPQ
        NH = nAlignment;
        if (nAlignment == 0) return;
        else if (nAlignment == 1) MAPQ = "60";
        else MAPQ = "1";
    }

    void extractFlagInfo() {
        paired = (flag & 1) != 0;
        forward = (flag & 16) == 0;
        if ((flag & 256) == 0) { // change all read to secondary alignment
            flag += 256;
        }
        primaryAlignment = false;
        mapped = (flag & 4) == 0;
        if (flag & 128) {
            pairSegment = 1;
        } else {
            pairSegment = 0; // it could be the first pair segment or it is single read.
        }
        concordant = (flag & 2) != 0;

        if (!mapped) {
            repeat = false;
        } else {
            if (chromosomeName.toZBuf() == "=") {
                pairToRepeat = true;
            }
        }
    }

    void praseCigar() {
        // input the Cigar string from SAM information and split it into each segment.
        // this function run faster than regex.

        cigarSegments.clear();
        cigarLength = 0;
        string cigar_string = cigarString.toZBuf();
        int previousLocation = 0;
        int lenLength = 0;
        for (int i = 0; i < cigar_string.size(); i++) {
            if (isalpha(cigar_string[i])){
                int len = stoi(cigar_string.substr(previousLocation, lenLength));
                cigarLength += len;
                cigarSegments.push_back(Cigar(len, cigar_string[i]));
                previousLocation = i + 1;
                lenLength = 0;
            } else {
                lenLength++;
            }
        }
    }

    void updateRA(char& read, char& ref) {
        // update RA tag
        RA_Array[asc2dna[ref]][asc2dna[read]]++;
    }

    void updateMP(int readLocation, char read, long long int refLocation, char ref) {
        char buf[1024];
        if (!MP.empty()) {
            MP.append(",");
        }
        itoa10<int>(asc2dna[ref] * 5 + asc2dna[read], buf);
        MP.append(buf);
        MP.append(":");
        itoa10<int>(readLocation + 1, buf);
        MP.append(buf);
        MP.append(":");
        itoa10<int>(refLocation + 1, buf);
        MP.append(buf);
    }

    void clearExtraTags() {
        TC = 0;
        MP.clear();
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                RA_Array[i][j] = 0;
            }
        }
    }

    void makeYZ(BTString &YZ_string) {
        if (TC == 0) {
            // did not find any conversion
            if (forward) {
                YZ_string = "OT";
            } else {
                YZ_string = "OB";
            }
        } else {
            if (forward && (conversionCount[0]<conversionCount[1])) {
                YZ_string = "OT";
            } else if (forward && (conversionCount[0]>=conversionCount[1])) {
                YZ_string = "CTOB";
            } else if (!forward && (conversionCount[0]<conversionCount[1])) {
                YZ_string = "CTOT";
            } else {
                YZ_string = "OB";
            }
        }
    }

    bool constructRepeatMD(BitPairReference* bitReference, MappingPositions &existPositions, bool &multipleAlignmentDetected) {

        praseCigar();

        ht2_error_t err = ht2_repeat_expand((TLAcycle == 0 || TLAcycle == 3) ? repeatHandles[0] : repeatHandles[1],
                                            chromosomeName.toZBuf(),
                                            location - 1,
                                            readSequence.length(),
                                            &repeatResult);

        BTString chromosomeRepeat;
        long long int locationRepeat;
        for (int i = 0; i < repeatResult->count; i++) {
            clearExtraTags();

            struct ht2_position *pos = &repeatResult->positions[i];
            chromosomeRepeat = refNameMap->names[pos->chr_id];
            for (int j = 0; j < chromosomeRepeat.length(); j++) {
                if (chromosomeRepeat[j] == ' ') {
                    chromosomeRepeat.trimEnd(chromosomeRepeat.length() - j);
                    break;
                }
            }
            locationRepeat = (pos->pos) + 1;
            bool genomeForward = pos->direction == 0;
            if (!genomeForward) {
                //cerr << readName << endl;
                continue;
            }

            /*if (!existPositions.append(chromosomeRepeat, locationRepeat)){
                continue;
            }*/

            ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
            SStringExpandable<char> raw_refbuf;
            raw_refbuf.resize(cigarLength + 16);
            raw_refbuf.clear();
            int off = bitReference->getStretch(
                    reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
                    (size_t)pos->chr_id,
                    (size_t)max<int>(locationRepeat-1, 0),
                    (size_t)cigarLength ASSERT_ONLY(, destU32));
            char* refSeq = raw_refbuf.wbuf() + off;

            BTString refSequence;
            refSequence.resize(cigarLength);
            for (int j = 0; j < cigarLength; j++) {
                //refSequence[j] = intToBase[*(refSeq + j)];
                refSequence.set(intToBase[*(refSeq + j)], j);
            }

            int repeatPositionsIndex;
            if (repeatPositions.sequenceExist(refSequence, repeatPositionsIndex)) {
                repeatPositions.appendRepeat(chromosomeRepeat, locationRepeat, repeatPositionsIndex);
                if (repeatPositions.positions[repeatPositionsIndex].AS == bestAS) {
                    nBestRepeat++;
                }
                continue;
            }

            BTString newMD;
            int nIncorrectMatch = 0;
            BTString repeatYZ;
            if (!constructRepeatMD(refSequence, newMD, nIncorrectMatch, repeatYZ)) {
                continue;
            }

            int newXM = XM + nIncorrectMatch;
            int newNM = NM + nIncorrectMatch;
            int newAS = AS - 6*nIncorrectMatch;

            repeatPositions.appendRepeat(locationRepeat, chromosomeRepeat, refSequence,newAS, newMD, newXM, newNM, TC, MP, RA_Array, repeatYZ);
            if (!paired) {
                if (newAS > bestAS) {
                    bestAS = newAS;
                    nBestRepeat = 1;
                } else if (newAS == bestAS) {
                    nBestRepeat++;
                }
                if (bestAS == 0 && nBestRepeat > 1 && uniqueOutputOnly) {
                    multipleAlignmentDetected = true;
                    return false;
                }
            }
            if (nBestRepeat > repeatLimit || i > repeatLimit) {
                if (nBestRepeat == 0 && !paired) {
                    continue;
                } else {
                    return true;
                }
            }
        }
        if ((nBestRepeat == 0 && !paired) || (paired && repeatPositions.positions.size() == 0)) {
            return false;
        } else {
            return true;
        }
    }

    bool constructRepeatMD(BTString &refSeq, BTString &newMD_String, int &nIncorrectMatch, BTString &repeatYZ) {
        // construct new MD string, this function is for repeat read
        // return true, if the read is mapped to correct location.
        // return false, if the read is mapped to incorrect location.

        char buf[1024];

        conversionCount[0] = 0;
        conversionCount[1] = 0;

        int readPos = 0;
        long long int refPos = 0;
        int count = 0;
        int newXM = 0;

        char cigarSymbol;
        int cigarLen;
        for (int i = 0; i < cigarSegments.size(); i++) {
            cigarSymbol = cigarSegments[i].getLabel();
            cigarLen = cigarSegments[i].getLen();

            if (cigarSymbol == 'S') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'N') {
                refPos += cigarLen;
            } else if (cigarSymbol == 'M') {
                for (int j = 0; j < cigarLen; j++) {
                    char readChar = readSequence[readPos];
                    char refChar = refSeq[refPos];
                    if (readChar == refChar || readChar == 'N' || refChar == 'N') {
                        count++;
                    } else {// mismatch
                        // output matched count
                        if (count != 0) {
                            //std::stringstream ss;
                            //ss << count;
                            //newMD_String.append(ss.str().c_str());
                            itoa10<int>(count, buf);
                            newMD_String.append(buf);
                            count = 0;
                        }
                        // output mismatch
                        if (!newMD_String.empty() && isalpha(newMD_String[newMD_String.length()-1])) {
                            char test = newMD_String[newMD_String.length()-1];
                            newMD_String.append('0');
                        }
                        if ((readChar == convertedTo) && (refChar == convertedFrom)) {
                            conversionCount[0]++;
                        } else if ((readChar == convertedToComplement) && (refChar == convertedFromComplement)) {
                            conversionCount[1]++;
                        } else {
                            // real mismatch
                            newXM++;
                        }
                        //newMD_String += refChar;
                        newMD_String.append(refChar);
                        updateMP(refPos, readChar, refPos, refChar);
                    }
                    if (newXM + (conversionCount[0] >= conversionCount[1] ? conversionCount[1]:conversionCount[0])> readSequence.length()/25) {
                        return false;
                    }
                    updateRA(readChar, refChar);
                    readPos++;
                    refPos++;
                }
            } else if (cigarSymbol == 'I') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'D') {
                //newMD_String += '^';
                newMD_String.append('^');
                for (int j = 0; j < cigarLen; j++) {
                    //newMD_String += refSeq[refPos];
                    newMD_String.append(refSeq[refPos]);
                    refPos++;
                }
            }
        }

        if (count != 0) {
            //newMD_String += to_string(count);
            //std::stringstream ss;
            //ss << count;
            //newMD_String.append(ss.str().c_str());
            itoa10<int>(count, buf);
            newMD_String.append(buf);
        }
        if (isalpha(newMD_String[0])) {
            // = '0' + newMD_String;
            newMD_String.insert('0', 0);
        }
        if (isalpha(newMD_String[newMD_String.length()-1])) {
            //newMD_String += '0';
            newMD_String.append('0');
        }

        newXM -= XM;
        nIncorrectMatch = (conversionCount[0] >= conversionCount[1]) ? conversionCount[1] : conversionCount[0];
        TC = (conversionCount[0] >= conversionCount[1]) ? conversionCount[0] : conversionCount[1];
        nIncorrectMatch += newXM;

        makeYZ(repeatYZ);
        return true;
    }

    bool constructMD(BitPairReference* bitReference) {
        // this funtction is for non-repeat alignment.
        // return true, if the read is mapped to correct location.
        // return false, if the alignment has too many mismatch.
        if (!mapped) {
            return true;
        }
        char buf[1024];
        praseCigar();
        //BTString newMD_String;
        MD.clear();
        //locationPointer = reference.getPointer(chromosomeName.toZBuf(), location-1);

        ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
        SStringExpandable<char> raw_refbuf;
        raw_refbuf.resize(cigarLength + 16);
        raw_refbuf.clear();
        int off = bitReference->getStretch(
                reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
                (size_t)chromosomeIndex,
                (size_t)max<int>(location-1, 0),
                (size_t)cigarLength ASSERT_ONLY(, destU32));
        char* refSeq = raw_refbuf.wbuf() + off;

        int readPos = 0;
        long long int refPos = 0;
        int count = 0;
        int newXM = 0;

        char cigarSymbol;
        int cigarLen;
        for (int i = 0; i < cigarSegments.size(); i++) {
            cigarSymbol = cigarSegments[i].getLabel();
            cigarLen = cigarSegments[i].getLen();

            if (cigarSymbol == 'S') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'N') {
                refPos += cigarLen;
            } else if (cigarSymbol == 'M') {
                for (int j = 0; j < cigarLen; j++) {
                    char readChar = readSequence[readPos];
                    char refChar = intToBase[*(refSeq + refPos)];
                    if (readChar == refChar || readChar == 'N' || refChar == 'N') {
                        count++;
                    } else {// mismatch
                        // output matched count
                        if (count != 0) {
                            //newMD_String += to_string(count);
                            //std::stringstream ss;
                            //ss << count;
                            //MD.append(ss.str().c_str());
                            itoa10<int>(count, buf);
                            MD.append(buf);
                            count = 0;
                        }
                        // output mismatch
                        if (!MD.empty() && isalpha(isalpha(MD[MD.length()-1]))) {
                            //newMD_String += '0';
                            MD.append('0');
                        }
                        if ((readChar == convertedTo) && (refChar == convertedFrom)) {
                            conversionCount[0]++;
                        } else if ((readChar == convertedToComplement) && (refChar == convertedFromComplement)) {
                            conversionCount[1]++;
                        } else {
                            // real mismatch
                            newXM++;
                        }
                        //newMD_String += refChar;
                        MD.append(refChar);
                        updateMP(refPos, readChar, refPos, refChar);
                    }
                    if (newXM + (conversionCount[0] >= conversionCount[1] ? conversionCount[1]:conversionCount[0])> readSequence.length()/25) {
                        return false;
                    }
                    updateRA(readChar, refChar);
                    readPos++;
                    refPos++;
                }
            } else if (cigarSymbol == 'I') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'D') {
                if (count != 0) {
                    //std::stringstream ss;
                    //ss << count;
                    //MD.append(ss.str().c_str());
                    itoa10<int>(count, buf);
                    MD.append(buf);
                    count = 0;
                }
                //newMD_String += '^';
                MD.append('^');
                for (int j = 0; j < cigarLen; j++) {
                    //newMD_String += *(refSeq + refPos);
                    MD.append(*(refSeq + refPos));
                    refPos++;
                }
            }
        }

        if (count != 0) {
            //newMD_String += to_string(count);
            //std::stringstream ss;
            //ss << count;
            //MD.append(ss.str().c_str());
            itoa10<int>(count, buf);
            MD.append(buf);
        }
        if (isalpha(MD[0])) {
            //newMD_String = '0' + newMD_String;
            MD.insert('0', 0);
        }
        if (isalpha(MD[MD.length()-1])) {
            //newMD_String += '0';
            MD.append('0');
        }

        int extraIncorrectMatch = 0;
        if (conversionCount[0] >= conversionCount[1]) {
            extraIncorrectMatch = conversionCount[1];
            TC = conversionCount[0];
        } else {
            extraIncorrectMatch = conversionCount[0];
            TC = conversionCount[1];
        }

        newXM -= XM; // find new mismatch;
        extraIncorrectMatch += newXM;

        makeYZ(YZ);
        NM += extraIncorrectMatch;
        XM += extraIncorrectMatch;
        AS = AS - 6*extraIncorrectMatch;
        //MD = newMD_String;
        return true;
    }

    void outputFlags(BTString& o) {
        // this function is for non-repeat or repeat-non-expand output.
        // for repeat-expand output, use outputRepeatFlags function.
        char buf[1024];
        if (mapped) {
            // AS
            o.append("AS:i:");
            itoa10<int>(AS, buf);
            o.append(buf);
            o.append('\t');
            // NH
            o.append("NH:i:");
            itoa10<int>(NH, buf);
            o.append(buf);
            o.append('\t');
            // XM
            o.append("XM:i:");
            itoa10<int>(XM, buf);
            o.append(buf);
            o.append('\t');
            // NM
            o.append("NM:i:");
            itoa10<int>(NM, buf);
            o.append(buf);
            o.append('\t');
            // MD
            o.append("MD:Z:");
            o.append(MD.toZBuf());
            o.append('\t');
            if (paired) {
                // YS
                o.append("YS:i:");
                itoa10<int>(YS, buf);
                o.append(buf);
                o.append('\t');
            }
        }
        o.append(unChangedTags.toZBuf());

        if (mapped && !repeat) {
            o.append('\t');
            // YZ
            o.append("YZ:Z:");
            o.append(YZ.toZBuf());
            o.append('\t');
            // TC
            o.append("TC:i:");
            itoa10<int>(TC, buf);
            o.append(buf);
            o.append('\t');
            // RA
            o.append("RA:i:");
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    itoa10<int>(RA_Array[i][j], buf);
                    o.append(buf);
                    if (!(i == 4 && j == 4)) {
                        o.append(",");
                    }
                }
            }
            o.append('\t');
            // MP
            if (!MP.empty()) {
                o.append("MP:Z:");
                o.append(MP.toZBuf());
            }
        }
    }

    void outputRepeatFlags(BTString& o, int &inputAS, int &inputXM, int &inputNM, BTString &inputMD,
            int &inputTC, int inputRA[5][5], BTString &inputMP, BTString &repeatYZ) {
        // this function is for repeat-expand alignment output.
        char buf[1024];
        // AS
        o.append("AS:i:");
        itoa10<int>(inputAS, buf);
        o.append(buf);
        o.append('\t');
        // NH
        o.append("NH:i:");
        itoa10<int>(NH, buf);
        o.append(buf);
        o.append('\t');
        // XM
        o.append("XM:i:");
        itoa10<int>(inputXM, buf);
        o.append(buf);
        o.append('\t');
        // NM
        o.append("NM:i:");
        itoa10<int>(inputNM, buf);
        o.append(buf);
        o.append('\t');
        // MD
        o.append("MD:Z:");
        o.append(inputMD.toZBuf());
        o.append('\t');
        // YS
        if (paired) {
            o.append("YS:i:");
            itoa10<int>(YS, buf);
            o.append(buf);
            o.append('\t');
        }
        // unchanged Tags
        o.append(unChangedTags.toZBuf());
        o.append('\t');
        // YZ
        o.append("YZ:i:");
        o.append(repeatYZ.toZBuf());
        o.append('\t');
        // TC
        o.append("TC:i:");
        itoa10<int>(inputTC, buf);
        o.append(buf);
        o.append('\t');
        // RA
        o.append("RA:i:");
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                itoa10<int>(inputRA[i][j], buf);
                o.append(buf);
                if (!(i == 4 && j == 4)) {
                    o.append(",");
                }
            }
        }
        o.append('\t');
        // MP
        if (!inputMP.empty()) {
            o.append("MP:Z:");
            o.append(inputMP.toZBuf());
        }
    }

    void outputSingleAlginemnt(BTString& o) {
        // this function is for non-repeat or repeat-non-expand alignment output.
        char buf[1024];
        // readName
        o.append(readName.toZBuf());
        o.append('\t');
        // flag, if it is primary alignment, -256
        itoa10<int>(flag-primaryAlignment*256, buf);
        o.append(buf);
        o.append('\t');
        // chromosome
        o.append(chromosomeName.toZBuf());
        o.append('\t');
        // location
        itoa10<int>(location, buf);
        o.append(buf);
        o.append('\t');
        //MAPQ
        o.append(MAPQ.toZBuf());
        o.append('\t');
        // cigar
        o.append(cigarString.toZBuf());
        o.append('\t');
        // pair to chromosome
        //if (concordant || (!oppositePairAddresses.empty() && !oppositePairAddresses[0]->mapped) || pairToRepeat || !mapped) {
        if (mapped && (!oppositePairAddresses.empty() && !oppositePairAddresses[0]->mapped)){
            o.append("=");
            o.append('\t');
        } else {
            o.append("0");
            o.append('\t');
        }

        //} else {
        //    o.append(pairToChromosome.toZBuf());
        //    o.append('\t');
        //}
        // pair to location
        if (concordant) {
            itoa10<int>(pairToLocation, buf);
            o.append(buf);
            o.append('\t');
        } else {
            itoa10<int>(pairToLocation, buf);
            o.append(buf);
            o.append('\t');
        }
        // pairing distance
        if (paired) {
            itoa10<int>(pairToLocation - location, buf);
            o.append(buf);
            o.append('\t');
        } else {
            o.append('0');
            o.append('\t');
        }
        // read sequence
        o.append(readSequence.toZBuf());
        o.append('\t');
        // read quality
        o.append(readQuality.toZBuf());
        o.append('\t');
        // tags
        outputFlags(o);
        o.append('\n');
    }

    void outputRepeatOppositePair(BTString& o, Alignment* oppositeAlignment, BTString &chromosome, long long int &location, long long int &oppoLocation, bool primary) {
        char buf[1024];

        for (int i = 0; i < oppositeAlignment->repeatPositions.positions.size(); i++) {
            MappingPosition *currentPosition = &oppositeAlignment->repeatPositions.positions[i];
            for (int j = 0; j < currentPosition->repeatPositions.size(); j++) {
                if (currentPosition->repeatPositions[j].location == location && currentPosition->repeatPositions[j].chromosome == chromosome) {
                    // readName
                    o.append(oppositeAlignment->readName.toZBuf());
                    o.append('\t');
                    // flag
                    itoa10<int>(oppositeAlignment->flag-(256*primary), buf);
                    o.append(buf);
                    o.append('\t');
                    // chromosome
                    o.append(chromosome.toZBuf());
                    o.append('\t');
                    // location
                    itoa10<int>(location, buf);
                    o.append(buf);
                    o.append('\t');
                    //MAPQ
                    o.append(oppositeAlignment->MAPQ.toZBuf());
                    o.append('\t');
                    // cigar
                    o.append(oppositeAlignment->cigarString.toZBuf());
                    o.append('\t');
                    // pair to chromosome
                    o.append("=");
                    o.append('\t');
                    // pair to location
                    itoa10<int>(oppoLocation, buf);
                    o.append(buf);
                    o.append('\t');
                    // pairing distance
                    itoa10<int>(oppoLocation-currentPosition->repeatPositions[j].location, buf);
                    o.append(buf);
                    o.append('\t');
                    // read sequence
                    o.append(oppositeAlignment->readSequence.toZBuf());
                    o.append('\t');
                    // read quality
                    o.append(oppositeAlignment->readQuality.toZBuf());
                    o.append('\t');
                    // tags
                    outputRepeatFlags(o, currentPosition->AS, currentPosition->XM, currentPosition->NM, currentPosition->MD,
                                      currentPosition->TC, currentPosition->RA_Array, currentPosition->MP, currentPosition->repeatPositions[j].YZ);
                    o.append('\n');
                    return;
                }
            }
        }
    }

    int outputRepeatAlignment(BTString& o) {
        char buf[1024];
        int nOutput = 0;
        for (int i = 0; i < repeatPositions.positions.size(); i++) {
            MappingPosition *currentPosition = &repeatPositions.positions[i];
            for (int j = 0; j < currentPosition->repeatPositions.size(); j++) {
                if (!paired && (currentPosition->AS == bestAS) ||
                    (paired && ((currentPosition->repeatPositions[j].pairScore == pairScore) || !oppositePairAddresses[0]->mapped))) {
                    if (!oppositePairAddresses.empty() && paired) {
                        nOutput++;
                        if (pairSegment == 0) {
                            // readName
                            o.append(readName.toZBuf());
                            o.append('\t');
                            // flag
                            itoa10<int>(flag-256*(nOutput==1 && primaryAlignment), buf);
                            o.append(buf);
                            o.append('\t');
                            // chromosome
                            o.append(currentPosition->repeatPositions[j].chromosome.toZBuf());
                            o.append('\t');
                            // location
                            itoa10<int>(currentPosition->repeatPositions[j].location, buf);
                            o.append(buf);
                            o.append('\t');
                            //MAPQ
                            o.append(MAPQ.toZBuf());
                            o.append('\t');
                            // cigar
                            o.append(cigarString.toZBuf());
                            o.append('\t');
                            // pair to chromosome
                            o.append('=');
                            o.append('\t');
                            // pair to location
                            if (paired) {
                                itoa10<int>(currentPosition->repeatPositions[j].pairToLocation, buf);
                                o.append(buf);
                                o.append('\t');
                            } else {
                                itoa10<int>(pairToLocation, buf);
                                o.append(buf);
                                o.append('\t');
                            }
                            // pairing distance
                            if (paired) {
                                itoa10<int>(currentPosition->repeatPositions[j].pairToLocation-currentPosition->repeatPositions[j].location, buf);
                                o.append(buf);
                                o.append('\t');
                            } else {
                                itoa10<int>(pairingDistance, buf);
                                o.append(buf);
                                o.append('\t');
                            }
                            // read sequence
                            o.append(readSequence.toZBuf());
                            o.append('\t');
                            // read quality
                            o.append(readQuality.toZBuf());
                            o.append('\t');
                            // tags
                            outputRepeatFlags(o, currentPosition->AS, currentPosition->XM, currentPosition->NM, currentPosition->MD,
                                              currentPosition->TC, currentPosition->RA_Array, currentPosition->MP, currentPosition->repeatPositions[j].YZ);
                            o.append('\n');
                            for (int k = 0; k < oppositePairAddresses.size(); k++) {
                                if (oppositePairAddresses[k]->repeat) {
                                    outputRepeatOppositePair(o, oppositePairAddresses[k], currentPosition->repeatPositions[j].chromosome, currentPosition->repeatPositions[j].pairToLocation, currentPosition->repeatPositions[j].location, (nOutput==1 && primaryAlignment));
                                } else {
                                    if (!oppositePairAddresses[k]->mapped) {
                                        oppositePairAddresses[k]->chromosomeName = currentPosition->repeatPositions[j].chromosome.toZBuf();
                                        oppositePairAddresses[k]->location = currentPosition->repeatPositions[j].location;
                                        oppositePairAddresses[k]->pairToLocation = currentPosition->repeatPositions[j].location;
                                    }
                                    oppositePairAddresses[k]->output(o);
                                }
                            }
                        } else {
                            for (int k = 0; k < oppositePairAddresses.size(); k++) {
                                if (oppositePairAddresses[k]->repeat) {
                                    outputRepeatOppositePair(o, oppositePairAddresses[k], currentPosition->repeatPositions[j].chromosome, currentPosition->repeatPositions[j].pairToLocation, currentPosition->repeatPositions[j].location, (nOutput==1 && primaryAlignment));
                                } else {
                                    if (!oppositePairAddresses[k]->mapped) {
                                        oppositePairAddresses[k]->chromosomeName = currentPosition->repeatPositions[j].chromosome.toZBuf();
                                        oppositePairAddresses[k]->location = currentPosition->repeatPositions[j].location;
                                        oppositePairAddresses[k]->pairToLocation = currentPosition->repeatPositions[j].location;
                                    }
                                    oppositePairAddresses[k]->output(o);
                                }
                            }
                            // readName
                            o.append(readName.toZBuf());
                            o.append('\t');
                            // flag
                            itoa10<int>(flag-256*(nOutput==1 && primaryAlignment), buf);
                            o.append(buf);
                            o.append('\t');
                            // chromosome
                            o.append(currentPosition->repeatPositions[j].chromosome.toZBuf());
                            o.append('\t');
                            // location
                            itoa10<int>(currentPosition->repeatPositions[j].location, buf);
                            o.append(buf);
                            o.append('\t');
                            //MAPQ
                            o.append(MAPQ.toZBuf());
                            o.append('\t');
                            // cigar
                            o.append(cigarString.toZBuf());
                            o.append('\t');
                            // pair to chromosome
                            o.append('=');
                            o.append('\t');
                            // pair to location
                            if (paired) {
                                itoa10<int>(currentPosition->repeatPositions[j].pairToLocation, buf);
                                o.append(buf);
                                o.append('\t');
                            } else {
                                itoa10<int>(pairToLocation, buf);
                                o.append(buf);
                                o.append('\t');
                            }
                            // pairing distance
                            if (paired) {
                                itoa10<int>(currentPosition->repeatPositions[j].pairToLocation-currentPosition->repeatPositions[j].location, buf);
                                o.append(buf);
                                o.append('\t');
                            } else {
                                itoa10<int>(pairingDistance, buf);
                                o.append(buf);
                                o.append('\t');
                            }
                            // read sequence
                            o.append(readSequence.toZBuf());
                            o.append('\t');
                            // read quality
                            o.append(readQuality.toZBuf());
                            o.append('\t');
                            // tags
                            outputRepeatFlags(o, currentPosition->AS, currentPosition->XM, currentPosition->NM, currentPosition->MD,
                                              currentPosition->TC, currentPosition->RA_Array, currentPosition->MP, currentPosition->repeatPositions[j].YZ);
                            o.append('\n');
                        }

                    } else {
                        // readName
                        o.append(readName.toZBuf());
                        o.append('\t');
                        // flag
                        itoa10<int>(flag-256*(nOutput==1 && primaryAlignment), buf);
                        o.append(buf);
                        o.append('\t');
                        // chromosome
                        o.append(currentPosition->repeatPositions[j].chromosome.toZBuf());
                        o.append('\t');
                        // location
                        itoa10<int>(currentPosition->repeatPositions[j].location, buf);
                        o.append(buf);
                        o.append('\t');
                        //MAPQ
                        o.append(MAPQ.toZBuf());
                        o.append('\t');
                        // cigar
                        o.append(cigarString.toZBuf());
                        o.append('\t');
                        // pair to chromosome
                        o.append('*');
                        o.append('\t');
                        // pair to location
                        o.append('0');
                        o.append('\t');
                        // pairing distance
                        o.append('0');
                        o.append('\t');
                        // read sequence
                        o.append(readSequence.toZBuf());
                        o.append('\t');
                        // read quality
                        o.append(readQuality.toZBuf());
                        o.append('\t');
                        // tags
                        outputRepeatFlags(o, currentPosition->AS, currentPosition->XM, currentPosition->NM, currentPosition->MD,
                                          currentPosition->TC, currentPosition->RA_Array, currentPosition->MP, currentPosition->repeatPositions[j].YZ);
                        o.append('\n');
                    }
                }
            }
        }
        return nOutput;
    }

    int output(BTString& o) {
        int nOutput;
        if (paired) {
            if (oppositePairAddresses[0]->primaryAlignment) {
                primaryAlignment = true;
            } else {
                primaryAlignment = false;
            }
        }
        if (repeat) {
            if (pairSegment == 0 || (pairSegment == 1 && !pairToRepeat)) {
                nOutput = outputRepeatAlignment(o);
            } else {
                nOutput = 0;
            }
        } else {
            outputSingleAlginemnt(o);
            nOutput = 1;
        }
        return nOutput;
    }
};

class Alignments {
public:
    EList<Alignment*> alignments; // pool to store current alignment result.
    queue<Alignment*> freeAlignments; // free pointer pool for new alignment result. after output a alignment, return the pointer back to this pool.

    TReadId previousReadID;
    int bestAS; // for single-end read selection
    int bestPairScore; // for paired-end read selection;
    int nBestPair; // NH for paired-end output.
    bool concordantAlignmentExist; // for paired-end output, if there is concordant alignment, output the concordant one.
    MappingPositions existPositions; // check if same alignment result is exist. same alignment result could made by planA alignment and planB alignment.

    BTString readName[2]; // the read name could be different for segment 1 and segment 2.
    BTString readSequence[2]; // save the read sequence for output.
    BTString qualityScore[2]; // save the quality score for output.

    bool working;
    bool paired;
    int repeatPoolLimit = 10;
    bool multipleAligned = false;
    const int maxPairScore = 500000;

    BitPairReference* bitReference;
    bool DNA = false;
    int nRepeatAlignment = 0;

    void initialize() {
        bestAS = numeric_limits<int>::min();
        bestPairScore = numeric_limits<int>::min()+1;
        existPositions.initialize();
        working = false;
        nBestPair = 0;
        concordantAlignmentExist = false;
        paired = false;
        multipleAligned = false;
        nRepeatAlignment = 0;
        for (int i = 0; i < 2; i++) {
            readName[i].clear();
            readSequence[i].clear();
            qualityScore[i].clear();
        }
        while(!alignments.empty()){
            alignments[0]->initialize();
            freeAlignments.push(alignments[0]);
            alignments.erase(0);
        }
    }

    Alignments(BitPairReference* ref, bool inputDNA): bitReference(ref), DNA(inputDNA) {
        initialize();
        for (int i = 0; i < 10; i++) {
            Alignment *newAlignment = new Alignment();
            freeAlignments.push(newAlignment);
        }
    }
    ~Alignments() {
        while (!freeAlignments.empty()) {
            delete freeAlignments.front();
            freeAlignments.pop();
        }
        for (int i = 0; i < alignments.size(); i++) {
            delete alignments[i];
        }
    }

    void getFreeAlignmentPointer(Alignment* &newAlignment) {
        //Alignment* newAlignment;
        if (!freeAlignments.empty()) {
            newAlignment = freeAlignments.front();
            freeAlignments.pop();
        } else {
            newAlignment = new Alignment();
        }
    }

    bool alignmentExist(Alignment *newAlignment) {

        if (newAlignment->repeat) {
            return false;
        }

        BTString chromosome = newAlignment->chromosomeName;
        long long int location = newAlignment->location;

        for (int i = 0; i < alignments.size(); i++) {
            if (location == alignments[i]->location && chromosome == alignments[i]->chromosomeName) {
                return true;
            }
        }
        return false;
    }

    bool isConcordant(long long int &location1, bool &forward1, long long int &location2, bool &forward2) {
        if (location1 < location2) {
            if (forward1 == true && forward2 == false) {
                return true;
            }
        } else {
            if (forward1 == false && forward2 == true) {
                return true;
            }
        }
        return false;
    }

    bool isConcordant(Alignment* alignment1, Alignment* alignment2) {
        if (!alignment1->mapped || !alignment2->mapped) {
            return false;
        }
        if (alignment1->repeat || alignment2->repeat) {
            return true;
        }
        if (alignment1->location < alignment2->location) {
            if (alignment1->forward == true && alignment2->forward == false) {
                return true;
            }
        } else {
            if (alignment1->forward == false && alignment2->forward == true) {
                return true;
            }
        }
        return false;
    }


    int calculatePairScore(Alignment *alignment1, Alignment *alignment2, int &nPair) {
        // calculate the pairScore for a pair of alignment result. Output pair Score and number of pair (1 for non-repeat).
        // Do not update their pairScore.

        int pairScore = numeric_limits<int>::min();
        nPair = 0;
        if (alignment1->pairSegment == alignment2->pairSegment){
            // when 2 alignment results are from same pair segment, output the lowest score and number of pair equal zero.
            pairScore = numeric_limits<int>::min();
        } else if (!alignment1->mapped && !alignment2->mapped) {
            // both ummaped.
            pairScore = numeric_limits<int>::min()/2 - 1;
        } else if (!alignment1->mapped || !alignment2->mapped) {
            // one of the segment unmapped.
            pairScore = numeric_limits<int>::min()/2;
            if (alignment1->repeat) {
                for (int ii = 0; ii < alignment1->repeatPositions.positions.size(); ii++) {
                    nPair += alignment1->repeatPositions.positions[ii].repeatPositions.size();
                }
            } else if (alignment2->repeat) {
                for (int ii = 0; ii < alignment2->repeatPositions.positions.size(); ii++) {
                    nPair += alignment2->repeatPositions.positions[ii].repeatPositions.size();
                }
            } else {
                nPair = 1;
            }
        } else if ((alignment1->chromosomeName != alignment2->chromosomeName && ((alignment1->chromosomeName != alignment2->pairToChromosome) && (alignment1->pairToChromosome != alignment2->chromosomeName))) ||
                    ((alignment1->location != alignment2->pairToLocation) && (alignment1->pairToLocation != alignment2->location))){
            // both mapped, but they do not pair to each other.
            pairScore = numeric_limits<int>::min()/2 + 2;
            nPair = 1;
        } else if ((!alignment1->repeat && !alignment2->repeat)){
            // both mapped and (both non-repeat or not expand repeat)
            if (DNA) {
                pairScore = calculatePairScore_DNA(alignment1->location,
                                                   alignment1->AS,
                                                   alignment1->forward,
                                                   alignment2->location,
                                                   alignment2->AS,
                                                   alignment2->forward);
            } else {
                pairScore = calculatePairScore_RNA(alignment1->location,
                                                   alignment1->XM,
                                                   alignment1->forward,
                                                   alignment2->location,
                                                   alignment2->XM,
                                                   alignment2->forward);
            }

            if (pairScore > alignment1->pairScore) {
                alignment1->pairToLocation = alignment2->location;
            }
            if (pairScore > alignment2->pairScore) {
                alignment2->pairToLocation = alignment1->location;
            }
            nPair = 1;
        } else if (alignment1->repeat && !alignment2->repeat) {
            // 1st alignment is repeat, 2nd alignment is not repeat.
            int tmpPairScore;
            int AS2 = alignment2->AS;
            int XM2 = alignment2->XM;
            long long int location2 = alignment2->location;
            bool forward1 = alignment1->forward;
            bool forward2 = alignment2->forward;
            BTString chrmosome2 = alignment2->chromosomeName;

            for (int i = 0; i < alignment1->repeatPositions.positions.size(); i++) {
                MappingPosition *repeatPositions = &alignment1->repeatPositions.positions[i];
                for (int j = 0; j < repeatPositions->repeatPositions.size(); j++) {
                    RepeatPosition *repeatPosition = &repeatPositions->repeatPositions[j];
                    if(repeatPosition->chromosome != chrmosome2) {
                        continue;
                    }
                    if (DNA) {
                        tmpPairScore = calculatePairScore_DNA(repeatPosition->location,
                                                              repeatPositions->AS,
                                                              forward1,
                                                              location2,
                                                              AS2,
                                                              forward2);
                    } else {
                        tmpPairScore = calculatePairScore_RNA(repeatPosition->location,
                                                              repeatPositions->XM,
                                                              forward1,
                                                              location2,
                                                              XM2,
                                                              forward2);
                    }

                    if (tmpPairScore >= alignment1->pairScore) {
                        repeatPositions->pairScore = tmpPairScore;
                        repeatPosition->pairScore = tmpPairScore;
                        repeatPosition->pairToLocation = location2;
                    }
                    if (tmpPairScore >= alignment2->pairScore) {
                        alignment2->pairToLocation = repeatPosition->location;
                    }
                    if (tmpPairScore > pairScore) {
                        pairScore = tmpPairScore;
                        nPair = 1;
                    } else if (tmpPairScore == pairScore) {
                        nPair++;
                    }
                }
            }
        } else if (!alignment1->repeat && alignment2->repeat) {
            // 1st alignment is non-repeat, 2nd alignment is repeat alignment.

            int tmpPairScore;
            int AS1 = alignment1->AS;
            int XM1 = alignment1->XM;
            long long int location1 = alignment1->location;
            BTString chromosome1 = alignment1->chromosomeName;
            bool forward1 = alignment1->forward;
            bool forward2 = alignment2->forward;

            for (int i = 0; i < alignment2->repeatPositions.positions.size(); i++) {
                MappingPosition *repeatPositions = &alignment2->repeatPositions.positions[i];
                for (int j = 0; j < repeatPositions->repeatPositions.size(); j++) {
                    RepeatPosition *repeatPosition = &repeatPositions->repeatPositions[j];
                    if(chromosome1 != repeatPosition->chromosome) {
                        continue;
                    }
                    if (DNA) {
                        tmpPairScore = calculatePairScore_DNA(repeatPosition->location,
                                                              repeatPositions->AS,
                                                              forward2,
                                                              location1,
                                                              AS1,
                                                              forward1);
                    } else {
                        tmpPairScore = calculatePairScore_RNA(repeatPosition->location,
                                                              repeatPositions->XM,
                                                              forward2,
                                                              location1,
                                                              XM1,
                                                              forward1);
                    }

                    if (tmpPairScore >= alignment1->pairScore) {
                        alignment1->pairToLocation = repeatPosition->location;
                    }
                    if (tmpPairScore >= alignment2->pairScore) {
                        repeatPositions->pairScore = tmpPairScore;
                        repeatPosition->pairScore = tmpPairScore;
                        repeatPosition->pairToLocation = location1;
                    }
                    if (tmpPairScore > pairScore) {
                        pairScore = tmpPairScore;
                        nPair = 1;
                    } else if (tmpPairScore == pairScore) {
                        nPair++;
                    }
                }
            }
        } else { // both repeat
            pairScore = numeric_limits<int>::min()/2 + 1;
            int tmpPairScore;

            bool forward1 = alignment1->forward;
            bool forward2 = alignment2->forward;

            for (int i1 = 0; i1 < alignment1->repeatPositions.positions.size(); i1++) {
                MappingPosition *repeatPositions1 = &alignment1->repeatPositions.positions[i1];
                for (int j1 = 0; j1 < repeatPositions1->repeatPositions.size(); j1++) {
                    RepeatPosition *repeatPosition1 = &repeatPositions1->repeatPositions[j1];

                    for (int i2 = 0; i2 < alignment2->repeatPositions.positions.size(); i2++) {
                        MappingPosition *repeatPositions2 = &alignment2->repeatPositions.positions[i2];
                        for (int j2 = 0; j2 < repeatPositions2->repeatPositions.size(); j2++) {
                            RepeatPosition *repeatPosition2 = &repeatPositions2->repeatPositions[j2];
                            if(repeatPosition1->chromosome != repeatPosition2->chromosome) {
                                continue;
                            }
                            if (DNA) {
                                tmpPairScore = calculatePairScore_DNA(repeatPosition1->location,
                                                                      repeatPositions1->AS,
                                                                      forward1,
                                                                      repeatPosition2->location,
                                                                      repeatPositions2->AS,
                                                                      forward2);
                            } else {
                                tmpPairScore = calculatePairScore_RNA(repeatPosition1->location,
                                                                      repeatPositions1->XM,
                                                                      forward1,
                                                                      repeatPosition2->location,
                                                                      repeatPositions2->XM,
                                                                      forward2);
                            }
                            if (tmpPairScore > pairScore) {
                                pairScore = tmpPairScore;
                                nPair = 1;
                            } else if (tmpPairScore == pairScore) {
                                if ((alignment1->pairSegment == 0 && repeatPosition1->pairScore == tmpPairScore) ||
                                        (alignment2->pairSegment == 0 && repeatPosition2->pairScore == tmpPairScore)){

                                } else {
                                    nPair++;
                                }
                            }
                            if (tmpPairScore >= repeatPositions1->pairScore) {
                                repeatPositions1->pairScore = tmpPairScore;
                                repeatPosition1->pairToLocation = repeatPosition2->location;
                                repeatPosition1->pairScore = tmpPairScore;
                            }
                            if (tmpPairScore >= repeatPositions2->pairScore) {
                                repeatPositions2->pairScore = tmpPairScore;
                                repeatPosition2->pairScore = tmpPairScore;
                                repeatPosition2->pairToLocation = repeatPosition1->location;
                            }
                        }
                    }
                }
            }
        }
        return pairScore;
    }

    int calculatePairScore_DNA(long long int &location1, int &AS1, bool &forward1, long long int &location2, int &AS2, bool &forward2) {
        // this is the basic function to calculate pair score.
        // if the distance between 2 alignment is more than 1000, we reduce the score by the distance/10.
        // if two alignment is concordant we add 500000 to make sure to select the concordant pair as best pair.
        int score = 100*AS1 + 100*AS2;
        int distance = abs(location1 - location2);
        if (distance > 1000) {
            score -= distance/100;
        }
        if (distance > 500000) {
            return numeric_limits<int>::min();
        }
        if (isConcordant(location1, forward1, location2, forward2)) {
            score += 500000;
        }
        return score;
    }

    int calculatePairScore_RNA(long long int &location1, int &XM1, bool &forward1, long long int &location2, int &XM2, bool &forward2) {
        // this is the basic function to calculate pair score.
        // if the distance between 2 alignment is more than 100,000, we reduce the score by the distance/10.
        // if two alignment is concordant we add 500,000 to make sure to select the concordant pair as best pair.
        int score = -100*XM1 + -100*XM2;
        int distance = abs(location1 - location2);
        if (distance > 100000) {
            score -= distance/1000;
        }
        if (distance > 500000) {
            return numeric_limits<int>::min();
        }
        if (isConcordant(location1, forward1, location2, forward2)) {
            score += 500000;
        }
        return score;
    }

    int calculateNumAlignment() {
        // this function is to calculate number of alignment (NH) for single-end output.
        int nAlignment = 0;

        for (int i = 0; i < alignments.size(); i++) {
            if (alignments[i]->repeat) {
                nAlignment += alignments[i]->nBestRepeat;
            } else {
                if (alignments[i]->mapped) {
                    nAlignment++;
                }
            }
        }

        return nAlignment;
    }

    int calculateNumBestPair() {
        // this is a backup function to calculate the number of best pair (NH). this can be use in output_paired function.
        int numBestPair = 0;
        for (int i = 0; i < alignments.size(); i++) {
            if ((alignments[i]->pairScore == bestPairScore) && (alignments[i]->concordant == concordantAlignmentExist)) {
                if (alignments[i]->pairSegment == 0) {
                    if (alignments[i]->repeat) {
                        MappingPositions repeatPositions = alignments[i]->repeatPositions;
                        for (int ii = 0; ii < repeatPositions.positions.size(); ii++) {
                            MappingPosition *currentPosition = &repeatPositions.positions[ii];
                            for (int jj = 0; jj < currentPosition->repeatPositions.size(); jj++) {
                                if (paired && (currentPosition->repeatPositions[jj].pairScore == bestPairScore)) {
                                    numBestPair++;
                                }
                            }
                        }
                    } else {
                        numBestPair++;
                    }
                }
            }
        }
        return numBestPair;
    }


    void addNewAlignment_single(Alignment *newAlignment) {
        // add new single-end alignment into pool.
        working = true;

        // save read name, sequence, and quality score for output.
        if (readName[0].empty()) {
            readName[0] = newAlignment->readName;
        }
        if (readSequence[0].empty()) {
            readSequence[0] = newAlignment->originalFw;
        }
        if (qualityScore[0].empty()) {
            qualityScore[0] = newAlignment->readQualityFw;
        }

        if (uniqueOutputOnly && multipleAligned) {
            newAlignment->initialize();
            freeAlignments.push(newAlignment);
            working = false;
            return;
        }

        // check if alignment result exist
        if(!existPositions.append(newAlignment)) {
            newAlignment->initialize();
            freeAlignments.push(newAlignment);
            working = false;
            return;
        }

        //newAlignment->loadRefInto(bitReference, refChromosomeNames);

        paired = newAlignment->paired;
        int newAlignmentAS;
        if (newAlignment->repeat) {
            if (!newAlignment->constructRepeatMD(bitReference, existPositions, multipleAligned)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
            newAlignmentAS = newAlignment->bestAS;
        } else {
            // check mismatch, update tags
            if (!newAlignment->constructMD(bitReference)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
            newAlignmentAS = newAlignment->AS;
        }

        if (newAlignmentAS > bestAS) {
            // if new alignment has best AS, clear alignments, update best AS, push new alignment to alignments.
            if (!alignments.empty()) {
                for (int i = alignments.size() - 1; i >= 0; i--) {
                    alignments[i]->initialize();
                    freeAlignments.push(alignments[i]);
                }
                alignments.clear();
            }
            bestAS = newAlignment->AS;
            alignments.push_back(newAlignment);
        } else if (newAlignmentAS == bestAS) {
            // if new alignment has same AS to best AS, add it to alignments.
            if (alignmentExist(newAlignment)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
            } else {
                alignments.push_back(newAlignment);
            }
        } else {
            // if new alignment has lower AS, return it to free alignment pool.
            newAlignment->initialize();
            freeAlignments.push(newAlignment);
        }
        if (bestAS == 0 && alignments.size() > 1) {
            multipleAligned = true;
        }
        working = false;
    }

    void addNewAlignment_paired(Alignment *newAlignment) {
        // add new paired-end alignment into pool.
        working = true;

        // save read name, sequence, quality score information for output.
        int pairSegment = newAlignment->pairSegment;
        if (readName[pairSegment].empty()) {
            readName[pairSegment] = newAlignment->readName;
        }
        if (readSequence[pairSegment].empty()) {
            readSequence[pairSegment] = newAlignment->originalFw;
        }
        if (qualityScore[pairSegment].empty()) {
            qualityScore[pairSegment] = newAlignment->readQualityFw;
        }

        // if there are too many repeat alignment result, ignore it. otherwise it will consume hugh amount of time.
        // this could happen when we use repeat index.
        if (newAlignment->repeat) {
            nRepeatAlignment++;
        }
        if (((nRepeatAlignment > repeatPoolLimit) || (uniqueOutputOnly && multipleAligned))) {
            newAlignment->initialize();
            freeAlignments.push(newAlignment);
            working = false;
            return;
        }

        // check if alignment result exist
        if(!existPositions.append(newAlignment)) {
            newAlignment->initialize();
            freeAlignments.push(newAlignment);
            working = false;
            return;
        }

        paired = newAlignment->paired;

        if (newAlignment->repeat) {
            if (!newAlignment->constructRepeatMD(bitReference, existPositions, multipleAligned)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
        } else {
            // check mismatch, update tags
            if (!newAlignment->constructMD(bitReference)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
        }

        int newPairScore = numeric_limits<int>::min();
        int nPair = 0;
        int pairTo = -1;

        if (alignments.empty() || newAlignment->pairSegment == 0) {
            // skip pairing process, directly push back
        } else{
            // find pairScore.
            if ((alignments.back()->location == newAlignment->pairToLocation) &&
                (alignments.back()->chromosomeName == newAlignment->chromosomeName) &&
                (alignments.back()->pairSegment != newAlignment->pairSegment)) {
                newPairScore = calculatePairScore(newAlignment, alignments.back(), nPair);
                if (nPair > 0) {
                    pairTo = alignments.size() - 1;
                }
            }
            /*for (int i = 0; i < alignments.size(); i++) {
                int tmpPairScore = calculatePairScore(newAlignment, alignments.back(), nPair);
                if (tmpPairScore > newPairScore) {
                    newPairScore = tmpPairScore;
                }
            }*/
        }
        // update pair score and nBest pair information in this class.
        if (newPairScore > bestPairScore) {
            bestPairScore = newPairScore;
            nBestPair = nPair;
        } else if (newPairScore == bestPairScore) {
            nBestPair += nPair;
        }
        if (bestPairScore == maxPairScore && nBestPair > 1) {
            multipleAligned = true;
        }
        // update pair score, pair address information for paired alignment.
        if (pairTo > -1) {
            bool concordant = isConcordant(newAlignment, alignments[pairTo]);
            if (newPairScore == bestPairScore) {
                concordantAlignmentExist = concordant;
            }
            if (newPairScore > newAlignment->pairScore) {
                newAlignment->pairScore = newPairScore;
                newAlignment->oppositePairAddresses.clear();
                newAlignment->oppositePairAddresses.push_back(alignments[pairTo]);
                newAlignment->pairToLocation = alignments[pairTo]->location;
                newAlignment->setConcordant(concordant);
            } else if (newPairScore == newAlignment->pairScore) {
                newAlignment->oppositePairAddresses.push_back(alignments[pairTo]);
                newAlignment->pairToLocation = alignments[pairTo]->location;
                newAlignment->setConcordant(concordant);
            }
            if (newPairScore > alignments[pairTo]->pairScore) {
                alignments[pairTo]->pairScore = newPairScore;
                alignments[pairTo]->oppositePairAddresses.clear();
                alignments[pairTo]->oppositePairAddresses.push_back(newAlignment);
                alignments[pairTo]->pairToLocation = newAlignment->location;
                alignments[pairTo]->setConcordant(concordant);
            } else if (newPairScore == alignments[pairTo]->pairScore) {
                alignments[pairTo]->oppositePairAddresses.push_back(newAlignment);
                alignments[pairTo]->pairToLocation = newAlignment->location;
                alignments[pairTo]->setConcordant(concordant);
            }
            alignments.push_back(newAlignment);
        } else {
            alignments.push_back(newAlignment);
        }
        working = false;
    }

    void outputUnAlignmentRead(BTString& o) {
        if (paired) {
            for (int i = 0; i < 2; i++) {
                string flag = (i == 0) ? "77" : "141";
                o.append(readName[i].toZBuf());
                o.append("\t");
                o.append(flag.c_str());
                o.append("\t*\t0\t0\t*\t*\t0\t0\t");
                o.append(readSequence[i].toZBuf());
                o.append("\t");
                o.append(qualityScore[i].toZBuf());
                o.append("\t");
                o.append("YT:Z:UP");
                o.append("\n");
            }
        } else {
            o.append(readName[0].toZBuf());
            o.append("\t4\t*\t0\t0\t*\t*\t0\t0\t");
            o.append(readSequence[0].toZBuf());
            o.append("\t");
            o.append(qualityScore[0].toZBuf());
            o.append("\t");
            o.append("YT:Z:UU");
            o.append("\n");
        }
    }

    void output_single(BTString& o, uint64_t &unAligned, uint64_t &nonRepeatAlignment, uint64_t &uniqueAlignment,
                uint64_t &multipleAlignment) {

        // find NH (nAlginment)
        int nAlignment = calculateNumAlignment();
        if (uniqueOutputOnly && multipleAligned) {
            nAlignment = 999;
        }
        // update statistics
        if (nAlignment == 0) {
            unAligned++;
        } else {
            nonRepeatAlignment++;
            if (nAlignment == 1) {
                uniqueAlignment++;
            } else {
                multipleAlignment++;
            }
        }
        // output
        bool primaryAlignment = true;
        int nOutput = 0;
        if (uniqueOutputOnly && nAlignment != 1) {
            // make a unalignment result and output it.
            //outputUnAlignmentRead(o);
            initialize();
            return;
        } else if (alignments.empty()) {
            // make a unalignment result and output it.
            outputUnAlignmentRead(o);
        } else {
            // update NH then output
            for (int i = 0; i < alignments.size(); i++) {
                alignments[i]->updateNH(nAlignment);
            }
            for (int i = 0; i < alignments.size(); i++) {
                if (primaryAlignment) {
                    alignments[i]->primaryAlignment = true;
                    primaryAlignment = false;
                }
                nOutput += alignments[i]->output(o);
            }
        }

        //assert (!uniqueOutputOnly && ((nOutput == nAlignment) || (nAlignment == 0 && nOutput == 1)));
        initialize();
    }

    void
    output_paired(BTString& o,
                  uint64_t &unConcordant,
                  uint64_t &uniqueConcordant,
                  uint64_t &multipleConcordant,
                  uint64_t &nonRepeatPairedAlignment,
                  uint64_t &uniqueDiscordant,
                  uint64_t &unAlignedPairRead,
                  uint64_t &alignedPairRead,
                  uint64_t &uniqueAlignedPairRead,
                  uint64_t &multipleAlignedPairRead) {

        //int testNumBestPair = calculateNumBestPair();

        // update statistics
        if (!concordantAlignmentExist) {
            unConcordant++;
            int nAlignment[2] = {0};
            if (nBestPair == 0) {
                unAlignedPairRead += 2;
            } else if (nBestPair == 1) {
                uniqueDiscordant++;
            } else {
                for (int i = 0; i < alignments.size(); i++) {
                    if (alignments[i]->pairScore == bestPairScore || alignments[i]->repeat) {
                        nAlignment[alignments[i]->pairSegment]++;
                    }
                }
                for (int i = 0; i < 2; i++) {
                    if (nAlignment[i] == 0) {
                        unAlignedPairRead++;
                    } else if (nAlignment[i] == 1) {
                        alignedPairRead++;
                        uniqueAlignedPairRead++;
                    } else {
                        alignedPairRead++;
                        multipleAlignedPairRead++;
                    }
                }
            }
        } else if (concordantAlignmentExist) {
            assert(nBestPair > 0);
            nonRepeatPairedAlignment++;
            if (nBestPair == 1) {
                uniqueConcordant++;
            } else {
                multipleConcordant++;
            }
        }

        // output
        int nOutput = 0;
        bool primaryAlignment = true;
        if (uniqueOutputOnly && (nBestPair != 1 || multipleAligned)) {
            // make a unalignment result and output it.
            //outputUnAlignmentRead(o);
            initialize();
            return;
        } else if (alignments.empty() || nBestPair == 0) {
            // make a unalignment result and output it.
            outputUnAlignmentRead(o);
        } else {
            // update NH then output
            for (int i = 0; i < alignments.size(); i++) {
                alignments[i]->updateNH(nBestPair);
            }
            for (int i = 0; i < alignments.size(); i++) {
                if ((alignments[i]->pairScore == bestPairScore) && (alignments[i]->concordant == concordantAlignmentExist)) {
                    if (primaryAlignment) {
                        alignments[i]->primaryAlignment = true;
                    }
                    // only output segment 0, then output it's pair (opposite pair).
                    if (alignments[i]->pairSegment == 0) {
                        if (alignments[i]->repeat) {
                            if (alignments[i]->primaryAlignment) {
                                alignments[i]->oppositePairAddresses[0]->primaryAlignment = true;
                            }
                            nOutput += alignments[i]->output(o);
                        } else {
                            for (int j = 0; j < alignments[i]->oppositePairAddresses.size(); j++) {
                                if (j > 0) {
                                    primaryAlignment = false;
                                    alignments[i]->primaryAlignment = false;
                                }
                                if (primaryAlignment) {
                                    alignments[i]->oppositePairAddresses[j]->primaryAlignment = true;
                                }
                                alignments[i]->oppositePairAddresses[j]->pairToLocation = alignments[i]->location;
                                alignments[i]->pairToLocation = alignments[i]->oppositePairAddresses[j]->location;
                                if (alignments[i]->repeat) {
                                    nOutput += alignments[i]->output(o);
                                } else if (alignments[i]->oppositePairAddresses[j]->repeat) {
                                    nOutput += alignments[i]->oppositePairAddresses[j]->output(o);
                                } else {
                                    if (!alignments[i]->mapped && alignments[i]->oppositePairAddresses[j]->mapped) {
                                        nOutput += alignments[i]->oppositePairAddresses[j]->output(o);
                                        alignments[i]->output(o);
                                    } else {
                                        nOutput += alignments[i]->output(o);
                                        alignments[i]->oppositePairAddresses[j]->output(o);
                                    }
                                }

                            }
                        }

                    }
                    if (alignments[i]->pairSegment == 1) {
                        primaryAlignment = false;
                    }
                }
            }
        }
        assert(nOutput == nBestPair);
        initialize();
    }


    void output(OutputQueue& oq_,
                size_t threadid_,
                uint64_t &unAligned,
                uint64_t &nonRepeatAlignment,
                uint64_t &uniqueAlignment,
                uint64_t &multipleAlignment,
                uint64_t &unConcordant,
                uint64_t &uniqueConcordant,
                uint64_t &multipleConcordant,
                uint64_t &nonRepeatPairedAlignment,
                uint64_t &uniqueDiscordant,
                uint64_t &unAlignedPairRead,
                uint64_t &alignedPairRead,
                uint64_t &uniqueAlignedPairRead,
                uint64_t &multipleAlignedPairRead) {

        while(working) {
            usleep(100);
        }
        BTString obuf_;
        OutputQueueMark qqm(oq_, obuf_, previousReadID, threadid_);
        if (paired) {
            output_paired(obuf_,
                          unConcordant,
                          uniqueConcordant,
                          multipleConcordant,
                          nonRepeatPairedAlignment,
                          uniqueDiscordant,
                          unAlignedPairRead,
                          alignedPairRead,
                          uniqueAlignedPairRead,
                          multipleAlignedPairRead);
        } else {
            output_single(obuf_,
                          unAligned,
                          nonRepeatAlignment,
                          uniqueAlignment,
                          multipleAlignment);
        }
    }
};

#endif //HISAT2_TLA_H