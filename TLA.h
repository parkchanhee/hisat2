//
// Created by Yun (Leo) Zhang on 1/10/20.
// This script is made for HISAT-3N (HISAT-TLA).
//

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
#include <unistd.h>
#include <queue>

class ReferenceGenome;
extern ReferenceGenome reference;
extern char convertedFrom;
extern char convertedTo;
extern char convertedFromComplement;
extern char convertedToComplement;
extern int incorrectTolerance;
extern bool expandRepeat;
extern vector<ht2_handle_t> repeatHandles;
extern struct ht2_index_getrefnames_result *refNameMap;
extern int repeatLimit;

using namespace std;

class Alignment;

class ChromosomeAddress{
    string chromosome;
    long long int startPosition;
public:
    ChromosomeAddress(string inputChromosome, long long int inputPosition):
            chromosome(inputChromosome), startPosition(inputPosition) {
    }

    const string& getChromosome() const {
        return chromosome;
    }

    long long int& getPosition() {
        return startPosition;
    }
};

class ReferenceGenome{
    // store reference genome in this class as a big character array for sequence and vector for chromosome start location

    char * refArray; // char array for reference genome storage.
    vector<ChromosomeAddress> chromosomeInfo; // vector to indicate start location (on refArray) for each chromosome.
    int N_chromosome = 0; // number of total chromosome. this is the right bondary for binary search.

    static bool sortChromosomeAddress(const ChromosomeAddress& first, const ChromosomeAddress& second) {
        // ascending sort rule by chromosome name
        return (first.getChromosome() < second.getChromosome());
    }

    void sortByChromosome() {
        // ascending sort by chromosome name
        sort(chromosomeInfo.begin(), chromosomeInfo.end(), sortChromosomeAddress);
    }

    long long int searchChromosome(string &chromosome, int left, int right) {
        // given chromosome name and search range, perform binary search to location chromosome in chromosomeInfo.

        if (left <= right) {
            int middle = (right + left) / 2;

            if (chromosomeInfo[middle].getChromosome() == chromosome) {
                return chromosomeInfo[middle].getPosition();
            }
            if (chromosomeInfo[middle].getChromosome() > chromosome) {
                return searchChromosome(chromosome, left, middle-1);
            }
            return searchChromosome(chromosome, middle+1, right);
        }

        throw "wrong chromosome name";
    }

    ifstream::pos_type getFileSize(string filename) {
        ifstream in(filename, std::ios::binary | std::ios::ate);
        return in.tellg();
    }

    vector<string> splitString(string &inputString, string delimiter) {
        vector<string> outputVector;
        size_t startPosition = 0;
        size_t endPosition = 0;
        string info;

        while ((endPosition = inputString.find(delimiter, startPosition)) != string::npos) {
            info = inputString.substr(startPosition, endPosition-startPosition);
            outputVector.push_back(info);
            startPosition = endPosition + delimiter.length();
        }
        outputVector.push_back(inputString.substr(startPosition));
        return outputVector;
    }
public:

    ReferenceGenome() {
    }

    ~ReferenceGenome() {
        free(refArray);
    }

    void load(string &referenceFileName) {
        // load genome reference sequence from disk to memory.
        // malloc a large memory with file size and write sequence into it.

        ifstream refFile;
        refFile.open(referenceFileName, ios_base::in);

        string chromosome;
        long long int n = 0; //this indicate current location in refArray.
        ifstream::pos_type filesize = getFileSize(referenceFileName);
        char* refArray_tmp = (char*) malloc(filesize);

        if (refFile.is_open()) {
            string line;

            while (refFile.good()) {
                getline(refFile, line);
                if (line.front() == '>') { // this line is chromosome name
                    chromosome = splitString(line, " ")[0].substr(1);
                    chromosomeInfo.push_back(ChromosomeAddress(chromosome, n));
                    N_chromosome = chromosomeInfo.size();
                } else {
                    for (int i = 0; i < line.size(); i++) {
                        *(refArray_tmp+n+i) = toupper(line[i]);
                    }
                    n += line.size();
                }
            }
        }

        refArray = (char*) realloc(refArray_tmp, n);
        sortByChromosome(); // sort my chromosome name for binary search.
        refFile.close();
    }

    char* getPointer(string chromosome, long long int location) {
        // by input chromosome name and the location on that chromosome, return the char pointer point to the location.
        return refArray + searchChromosome(chromosome, 0, N_chromosome-1) + location;
    }
};

class Cigar {
    int len;
    char label;
public:

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
    string chromosome;
    int pairScore;
    long long int pairToLocation;
    bool concordant;

    void initialize() {
        pairScore = numeric_limits<int>::min();
        pairToLocation = 0;
        concordant = false;
    }

    RepeatPosition(long long int &inputLocation, string &inputChromosome):location(inputLocation), chromosome(inputChromosome) {
        initialize();
    }
};

class MappingPosition {
public:
    long long int location;
    string chromosome;
    int AS;
    int pairScore;
    int pairSegment;
    long long int pairToLocation;
    bool concordant;

    // for repeat alignment only
    string MD;
    //int AS;
    int XM;
    int NM;
    int YS;
    int TC;
    BTString MP;
    int RA_Array[5][5] = {{0,},};
    string refSequence;
    vector<RepeatPosition> repeatPositions;

    MappingPosition (long long int inputLocation, string inputChromosome, int inputAS_score, int inputReadSegment=0, bool inputConcordant=false,
                     int inputPairScore=numeric_limits<int>::min(), long long int inputPairToLocation = -1){
        concordant = inputConcordant;
        location = inputLocation;
        chromosome = inputChromosome;
        AS = inputAS_score;
        pairScore = inputPairScore;
        pairSegment = inputReadSegment;
        pairToLocation = inputPairToLocation;
        pairScore = numeric_limits<int>::min();
    }

    MappingPosition (long long int inputLocation, string inputChromosome, string &inputRefSequence, int &inputAS, string &inputMD, int &inputXM, int &inputNM, int &inputTC, BTString &inputMP, int inputRA[5][5], int inputYS):
            location(inputLocation), chromosome(inputChromosome), refSequence(inputRefSequence), AS(inputAS), MD(inputMD), XM(inputXM), NM(inputNM), TC(inputTC), MP(inputMP), YS(inputYS){
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                RA_Array[i][j] = inputRA[i][j];
            }
        }
        pairScore = numeric_limits<int>::min();
        addRepeatPosition(inputLocation, inputChromosome);
    }

    void addRepeatPosition(long long int &inputLocation, string &inputChromosome) {
        repeatPositions.push_back(RepeatPosition(inputLocation, inputChromosome));
    }
};

class MappingPositions {
public:
    vector<MappingPosition> positions;
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

    bool sequenceExist (string& refSequence, int &index) {
        // return true if reference sequence is exist, else, return false.
        for (int i = 0; i < positions.size(); i++) {
            if (refSequence == positions[i].refSequence) {
                index = i;
                return true;
            }
        }
        return false;
    }

    bool positionExist (long long int &location, string &chromosome, int pairSegment, int &index) {
        // return true if position is exist, else, return false.
        for (int i = 0; i < positions.size(); i++) {
            if ((positions[i].pairSegment == pairSegment) &&
                (positions[i].location == location) &&
                (positions[i].chromosome == chromosome)) {
                index = i;
                return true;
            }
        }
        return false;
    }

    void appendRepeat(string &chromosome, long long int &location, int &index) {
        positions[index].addRepeatPosition(location, chromosome);
    }

    void appendRepeat(long long int &location, string &chromosome, string &refSequence, int &AS, string &MD, int &XM, int &NM, int &TC, BTString &MP, int RA[5][5], int YS=0) {
        positions.push_back(MappingPosition(location, chromosome, refSequence, AS, MD, XM, NM, TC, MP, RA, YS));
    }

    bool append(Alignment* newAlignment);

    bool append (string chromosome, long long int location, int AS, int pairSegment=0) {
        // return true if the position is not exist and will append to positions, else, false.
        int index;
        if (positionExist(location, chromosome, pairSegment, index)) {
            if (positions[index].AS == AS) {
                return false;
            } else {
                return true;
            }
        } else {
            positions.emplace_back(location, chromosome, AS);
            return true;
        }
    }
};

class Alignment {
public:
    BTString readName;
    int flag;
    BTString chromosomeName;
    long long int location;
    BTString MAPQ;
    BTString cigarString;
    BTString pairToChromosome;
    long long int pairToLocation;
    long long int pairingDistance;
    BTString readSequence;
    BTString readQuality;

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

    // unChanged tags
    BTString unChangedTags;

    // intermediate variable
    char* locationPointer; // point to the first char for mapping location in genome.
    bool planA;
    bool paired;
    bool forward;
    bool primaryAlignment; // this is for output flag adjustment.
    bool mapped;
    bool concordant;
    int pairSegment; // 0 for first segment, 1 for second segment.
    struct ht2_repeat_expand_result *repeatResult = nullptr;
    int pairScore; // to identify the better pair

    vector<Alignment*> oppositePairAddresses;
    // for repeatAlignment only
    bool repeat;
    bool pairToRepeat;
    MappingPositions repeatPositions; // only have chromosme and location informations
    int bestAS;
    int nBestRepeat;


    void initialize() {
        readName.clear();
        flag = -1;
        chromosomeName.clear();
        location = 0;
        MAPQ.clear();
        cigarString.clear();
        pairToChromosome.clear();
        pairToLocation = 0;
        pairingDistance = 0;
        readSequence.clear();
        readQuality.clear();

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

        if (repeatResult != nullptr) {
            free(repeatResult);
            repeatResult = nullptr;
        }
    }

    Alignment() {
        initialize();
    }

    ~Alignment() {
        if (repeatResult != nullptr) {
            free(repeatResult);
        }
    }

    void setConcordant(bool concordant_) {
        // change concordant status and flag
        concordant = concordant_;
        if ((flag&2) && !concordant) {
            flag -= 2;
        } else if (!(flag&2) && concordant) {
            flag += 2;
        }
    }
    void updateNH(int nAlignment) {
        if (!mapped) {
            return;
        }
        // update NH and MAPQ
        NH = nAlignment;
        if (nAlignment == 0) {
            return;
        } else if (nAlignment == 1) {
            MAPQ = "60";
        } else {
            MAPQ = "1";
        }
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

    int RA_Map(char &base) {
        // this is a helper function for RA tag construction.
        if (base == 'A') {
            return 0;
        }
        if (base == 'C') {
            return 1;
        }
        if (base == 'G') {
            return 2;
        }
        if (base == 'T') {
            return 3;
        }
        return 4;
    }

    vector<Cigar> getCigarSegments(string cigarString) {
        // input the Cigar string from SAM information and split it into each segment.
        // this function run faster than regex.
        vector<Cigar> outputCigar;
        int previousLocation = 0;
        int lenLength = 0;
        for (int i = 0; i < cigarString.size(); i++) {
            if (isalpha(cigarString[i])){
                outputCigar.push_back(Cigar( stoi(cigarString.substr(previousLocation, lenLength)),cigarString[i]));
                previousLocation = i + 1;
                lenLength = 0;
            } else {
                lenLength++;
            }
        }
        return outputCigar;
    }

    void updateRA(char& read, char& ref) {
        // update RA tag
        RA_Array[RA_Map(ref)][RA_Map(read)]++;
    }

    void updateMP(int readLocation, char read, long long int refLocation, char ref) {
        char buf[1024];
        if (!MP.empty()) {
            MP.append(",");
        }
        itoa10<int>(RA_Map(ref) * 5 + RA_Map(read), buf);
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

    bool constructRepeatMD(MappingPositions &existPositions) {

        string readSequence_string = readSequence.toZBuf();
        int sequenceLength = readSequence_string.size();

        ht2_error_t err = ht2_repeat_expand(planA ? repeatHandles[0] : repeatHandles[1],
                                            chromosomeName.toZBuf(),
                                            location - 1,
                                            sequenceLength,
                                            &repeatResult);

        string chromosomeRepeat;
        long long int locationRepeat;
        for (int i = 0; i < repeatResult->count; i++) {
            clearExtraTags();

            struct ht2_position *pos = &repeatResult->positions[i];
            chromosomeRepeat = refNameMap->names[pos->chr_id];
            if (chromosomeRepeat.find(' ', 0) != string::npos) {
                chromosomeRepeat = chromosomeRepeat.substr(0, chromosomeRepeat.find(' ', 0));
            }
            bool genomeForward = pos->direction == 0;
            if (!genomeForward) {
                continue;
            }
            locationRepeat = (pos->pos) + 1;

            if (!existPositions.append(chromosomeRepeat, locationRepeat, AS)){
                continue;
            }

            locationPointer = reference.getPointer(chromosomeRepeat, locationRepeat-1);
            string refSequence(locationPointer, sequenceLength);

            int repeatPositionsIndex;
            if (repeatPositions.sequenceExist(refSequence, repeatPositionsIndex)) {
                repeatPositions.appendRepeat(chromosomeRepeat, locationRepeat, repeatPositionsIndex);
                if (repeatPositions.positions[repeatPositionsIndex].AS == bestAS) {
                    nBestRepeat++;
                }
                continue;
            }

            string newMD;
            int nIncorrectMatch = 0;
            if (!constructRepeatMD(readSequence_string, newMD, nIncorrectMatch)) {
                continue;
            }

            int newXM = XM - nIncorrectMatch;
            int newNM = NM - nIncorrectMatch;
            int newAS = AS - 5*nIncorrectMatch;

            repeatPositions.appendRepeat(locationRepeat, chromosomeRepeat, refSequence,newAS, newMD, newXM, newNM, TC, MP, RA_Array);
            if (!paired) {
                if (newAS > bestAS) {
                    bestAS = newAS;
                    nBestRepeat = 1;
                } else if (newAS == bestAS) {
                    nBestRepeat++;
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
        if (nBestRepeat == 0 && !paired) {
            return false;
        } else {
            return true;
        }
    }

    bool constructRepeatMD(string &readSequence_string, string &newMD_String, int &nIncorrectMatch) {
        // construct new MD string, this function is for repeat read
        // return true, if the read is mapped to correct location.
        // return false, if the read is mapped to incorrect location.

        vector<Cigar> cigarSegments = getCigarSegments(cigarString.toZBuf());

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
                    char readChar = readSequence_string[readPos];
                    char refChar = *(locationPointer + refPos);
                    if (readChar == refChar || readChar == 'N' || refChar == 'N') {
                        count++;
                    } else {// mismatch
                        // output matched count
                        if (count != 0) {
                            newMD_String += to_string(count);
                            count = 0;
                        }
                        // output mismatch
                        if (!newMD_String.empty() && isalpha(newMD_String.back())) {
                            newMD_String += '0';
                        }
                        if (((readChar == convertedTo) && (refChar == convertedFrom) && planA) ||
                            ((readChar == convertedToComplement) && (refChar == convertedFromComplement) && !planA)) {
                            // converted base, output to MD, do not count as mismatch
                            newMD_String += refChar;
                            TC++;
                        } else {
                            // real mismatch
                            newMD_String += refChar;
                            newXM++;
                            if (newXM - XM > incorrectTolerance) {
                                return false;
                            }
                        }
                        updateMP(refPos, readChar, refPos, refChar);
                    }
                    updateRA(readChar, refChar);
                    readPos++;
                    refPos++;
                }
            } else if (cigarSymbol == 'I') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'D') {
                newMD_String += '^';
                for (int j = 0; j < cigarLen; j++) {
                    newMD_String += *(locationPointer + refPos);
                    refPos++;
                }
            }
        }

        if (count != 0) {
            newMD_String += to_string(count);
        }
        if (isalpha(newMD_String.front())) {
            newMD_String = '0' + newMD_String;
        }
        if (isalpha(newMD_String.back())) {
            newMD_String += '0';
        }

        return true;
    }

    bool constructMD() {
        // this funtction is for non-repeat alignment.
        // return true, if the read is mapped to correct location.
        // return false, if the alignment has too many mismatch.
        if (!mapped) {
            return true;
        }
        vector<Cigar> cigarSegments = getCigarSegments(cigarString.toZBuf());
        string newMD_String;
        locationPointer = reference.getPointer(chromosomeName.toZBuf(), location-1);

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
                    char refChar = *(locationPointer + refPos);
                    if (readChar == refChar || readChar == 'N' || refChar == 'N') {
                        count++;
                    } else {// mismatch
                        // output matched count
                        if (count != 0) {
                            newMD_String += to_string(count);
                            count = 0;
                        }
                        // output mismatch
                        if (!newMD_String.empty() && isalpha(newMD_String.back())) {
                            newMD_String += '0';
                        }
                        if (((readChar == convertedTo) && (refChar == convertedFrom) && planA) ||
                                ((readChar == convertedToComplement) && (refChar == convertedFromComplement) && !planA)) {
                            // converted base, output to MD, do not count as mismatch
                            newMD_String += refChar;
                            TC++;
                        } else {
                            // real mismatch
                            newMD_String += refChar;
                            newXM++;
                            if (newXM - XM > incorrectTolerance) {
                                return false;
                            }
                        }
                        updateMP(refPos, readChar, refPos, refChar);
                    }
                    updateRA(readChar, refChar);
                    readPos++;
                    refPos++;
                }
            } else if (cigarSymbol == 'I') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'D') {
                if (count != 0) {
                    newMD_String += to_string(count);
                    count = 0;
                }
                newMD_String += '^';
                for (int j = 0; j < cigarLen; j++) {
                    newMD_String += *(locationPointer + refPos);
                    refPos++;
                }
            }
        }

        if (count != 0) {
            newMD_String += to_string(count);
        }
        if (isalpha(newMD_String.front())) {
            newMD_String = '0' + newMD_String;
        }
        if (isalpha(newMD_String.back())) {
            newMD_String += '0';
        }

        int nIncorrectMatch = newXM - XM;
        NM += nIncorrectMatch;
        XM = newXM;
        AS = AS - 5*nIncorrectMatch;
        MD = newMD_String;
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

    void outputRepeatFlags(BTString& o, int &inputAS, int &inputXM, int &inputNM, string &inputMD,
            int &inputTC, int inputRA[5][5], BTString &inputMP) {
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
        o.append(inputMD.c_str());
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

    void outputRepeatOppositePair(BTString& o, Alignment* oppositeAlignment, string &chromosome, long long int &location, long long int &oppoLocation, bool primary) {
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
                    o.append(chromosome.c_str());
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
                                      currentPosition->TC, currentPosition->RA_Array, currentPosition->MP);
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
                            o.append(currentPosition->repeatPositions[j].chromosome.c_str());
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
                                              currentPosition->TC, currentPosition->RA_Array, currentPosition->MP);
                            o.append('\n');
                            for (int k = 0; k < oppositePairAddresses.size(); k++) {
                                if (oppositePairAddresses[k]->repeat) {
                                    outputRepeatOppositePair(o, oppositePairAddresses[k], currentPosition->repeatPositions[j].chromosome, currentPosition->repeatPositions[j].pairToLocation, currentPosition->repeatPositions[j].location, (nOutput==1 && primaryAlignment));
                                } else {
                                    if (!oppositePairAddresses[k]->mapped) {
                                        oppositePairAddresses[k]->chromosomeName = currentPosition->repeatPositions[j].chromosome.c_str();
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
                                        oppositePairAddresses[k]->chromosomeName = currentPosition->repeatPositions[j].chromosome.c_str();
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
                            o.append(currentPosition->repeatPositions[j].chromosome.c_str());
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
                                              currentPosition->TC, currentPosition->RA_Array, currentPosition->MP);
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
                        o.append(currentPosition->repeatPositions[j].chromosome.c_str());
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
                                          currentPosition->TC, currentPosition->RA_Array, currentPosition->MP);
                        o.append('\n');
                    }
                }
            }
        }
        return nOutput;
    }

    int output(BTString& o) {
        int nOutput;
        if (oppositePairAddresses[0]->primaryAlignment) {
            primaryAlignment = true;
        } else {
            primaryAlignment = false;
        }
        if (repeat && expandRepeat) {
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
    vector<Alignment*> alignments; // pool to store current alignment result.
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

    void initialize() {
        bestAS = numeric_limits<int>::min();
        bestPairScore = numeric_limits<int>::min()+1;
        existPositions.initialize();
        working = false;
        nBestPair = 0;
        concordantAlignmentExist = false;
        paired = false;
        for (int i = 0; i < 2; i++) {
            readName[i].clear();
            readSequence[i].clear();
            qualityScore[i].clear();
        }
        while(!alignments.empty()){
            alignments[0]->initialize();
            freeAlignments.push(alignments[0]);
            alignments.erase(alignments.begin());
        }
        while (freeAlignments.size() > 30) {
            delete freeAlignments.front();
            freeAlignments.pop();
        }
    }

    Alignments() {
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
        } else if (alignment1->planA != alignment2->planA) {
            // both mapped, but they belong to different plan.
            pairScore = numeric_limits<int>::min()/2;
            if (alignment1->chromosomeName == alignment2->chromosomeName) {
                pairScore += 1;
            }
            nPair = 1;
        }else if ((alignment1->chromosomeName != alignment2->chromosomeName && ((alignment1->chromosomeName != alignment2->pairToChromosome) && (alignment1->pairToChromosome != alignment2->chromosomeName))) ||
                    ((alignment1->location != alignment2->pairToLocation) && (alignment1->pairToLocation != alignment2->location))){
            // both mapped, but they do not pair to each other.
            pairScore = numeric_limits<int>::min()/2 + 2;
            nPair = 1;
        } else if ((!alignment1->repeat && !alignment2->repeat) || !expandRepeat){
            // both mapped and (both non-repeat or not expand repeat)
            pairScore = calculatePairScore(alignment1->location, alignment1->AS, alignment1->forward, alignment2->location, alignment2->AS, alignment2->forward);
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
            long long int location2 = alignment2->location;
            bool forward1 = alignment1->forward;
            bool forward2 = alignment2->forward;
            string chrmosome2 = alignment2->chromosomeName.toZBuf();

            for (int i = 0; i < alignment1->repeatPositions.positions.size(); i++) {
                MappingPosition *repeatPositions = &alignment1->repeatPositions.positions[i];
                for (int j = 0; j < repeatPositions->repeatPositions.size(); j++) {
                    RepeatPosition *repeatPosition = &repeatPositions->repeatPositions[j];
                    if(repeatPosition->chromosome != chrmosome2) {
                        continue;
                    }
                    tmpPairScore = calculatePairScore(repeatPosition->location, repeatPositions->AS, forward1, location2, AS2, forward2);
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
            long long int location1 = alignment1->location;
            string chromosome1 = alignment1->chromosomeName.toZBuf();
            bool forward1 = alignment1->forward;
            bool forward2 = alignment2->forward;

            for (int i = 0; i < alignment2->repeatPositions.positions.size(); i++) {
                MappingPosition *repeatPositions = &alignment2->repeatPositions.positions[i];
                for (int j = 0; j < repeatPositions->repeatPositions.size(); j++) {
                    RepeatPosition *repeatPosition = &repeatPositions->repeatPositions[j];
                    if(chromosome1 != repeatPosition->chromosome) {
                        continue;
                    }
                    tmpPairScore = calculatePairScore(repeatPosition->location, repeatPositions->AS, forward2, location1, AS1, forward1);
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
                            tmpPairScore = calculatePairScore(repeatPosition1->location, repeatPositions1->AS, forward1, repeatPosition2->location, repeatPositions2->AS, forward2);
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


    int calculatePairScore(long long int &location1, int &AS1, bool &forward1, long long int &location2, int &AS2, bool &forward2) {
        // this is the basic function to calculate pair score.
        // if the distance between 2 alignment is more than 1000, we reduce the score by the distance.
        // if two alignment is concordant we add 500000 to make sure to select the concordant pair as best pair.
        int score = 100*AS1 + 100*AS2;
        int distance = abs(location1 - location2);
        if (distance > 1000) {
            score -= distance;
        }
        if (distance > 500000) {
            return numeric_limits<int>::min();
        }
        if (isConcordant(location1, forward1, location2, forward2)) {
            score += 500000;
        }
        return score;
    }

    int calculateNumAlignment(bool& repeat) {
        // this function is to calculate number of alignment (NH) for single-end output.
        int nAlignment = 0;
        if (expandRepeat) {
            for (int i = 0; i < alignments.size(); i++) {
                if (alignments[i]->repeat) {
                    nAlignment += alignments[i]->nBestRepeat;
                    repeat = true;
                } else {
                    if (alignments[i]->mapped) {
                        nAlignment++;
                    }
                }
            }
        } else {
            for (int i = 0; i < alignments.size(); i++) {
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
                    if (alignments[i]->repeat && expandRepeat) {
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

    void addRepeatIndexTag(Alignment *newAlignment) {
        // add a extra tag if the read is aligned on repeat index. this tag will show which repeat index is used.
        // XR:Z:TC means the read aligned on TC conversion repeat index.
        string tag = "";
        if (newAlignment->planA) {
            tag += convertedFrom;
            tag += convertedTo;
        } else {
            tag += convertedFromComplement;
            tag += convertedToComplement;
        }
        if (!newAlignment->unChangedTags.empty()) {
            newAlignment->unChangedTags.append('\t');
        }
        newAlignment->unChangedTags.append("XR:Z:");
        newAlignment->unChangedTags.append(tag.c_str());
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

        paired = newAlignment->paired;
        int newAlignmentAS;
        if (newAlignment->repeat && expandRepeat) {
            if (!newAlignment->constructRepeatMD(existPositions)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
            newAlignmentAS = newAlignment->bestAS;
        } else {
            // check if alignment result exist
            if(!existPositions.append(newAlignment)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
            if (newAlignment->repeat) {
                addRepeatIndexTag(newAlignment);
            } else {
                // check mismatch, update tags
                if (!newAlignment->constructMD()) {
                    newAlignment->initialize();
                    freeAlignments.push(newAlignment);
                    working = false;
                    return;
                }
            }
            newAlignmentAS = newAlignment->AS;
        }

        if (newAlignmentAS > bestAS) {
            // if new alignment has best AS, clear alignments, update best AS, push new alignment to alignments.
            for (int i = 0; i < alignments.size(); i++) {
                alignments[i]->initialize();
                freeAlignments.push(alignments[i]);
            }
            alignments.clear();
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
        working = false;
    }

    void addNewAlignment_paired(Alignment *newAlignment) {
        // add new paired-end alignment into pool.
        working = true;
        // if there are too many alignment result, ignore it. otherwise it will consume hugh amount of time.
        // this could happen when we use repeat index.
        if (alignments.size() > 20) {
            newAlignment->initialize();
            freeAlignments.push(newAlignment);
            working = false;
            return;
        }
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
        paired = newAlignment->paired;

        if (newAlignment->repeat && expandRepeat) {
            // check if alignment result exist
            if(!existPositions.append(newAlignment)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
            if (!newAlignment->constructRepeatMD(existPositions)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
        } else {
            // check if alignment result exist
            if(!existPositions.append(newAlignment)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
                working = false;
                return;
            }
            if (newAlignment->repeat) {
                addRepeatIndexTag(newAlignment);
            } else {
                // check mismatch, update tags
                if (!newAlignment->constructMD()) {
                    newAlignment->initialize();
                    freeAlignments.push(newAlignment);
                    working = false;
                    return;
                }
            }
        }

        int bestNewPairScore = numeric_limits<int>::min();
        int tmp_nBestPair = 0;
        int pairTo = -1;

        if (alignments.empty()) {
            // jump pairing process, directly push back
        } else{
            for (int i = 0; i < alignments.size(); i++) {
                // find it's pair and pairScore.
                if ((alignments[i]->location != newAlignment->pairToLocation) || (alignments[i]->chromosomeName != newAlignment->chromosomeName)) {
                    continue;
                }
                if (alignments[i]->pairSegment == newAlignment->pairSegment) {
                    continue;
                }
                int nPair = 0;
                int newPairScore = calculatePairScore(newAlignment, alignments[i], nPair);
                if (newPairScore >= bestNewPairScore) {
                    bestNewPairScore = newPairScore;
                    tmp_nBestPair = nPair;
                    pairTo = i;
                }
            }
        }
        // update pair score and nBest pair information in this class.
        if (bestNewPairScore > bestPairScore) {
            bestPairScore = bestNewPairScore;
            nBestPair = tmp_nBestPair;
        } else if (bestNewPairScore == bestPairScore) {
            nBestPair += tmp_nBestPair;
        }
        // update pair score, pair address information for paired alignment.
        if (pairTo > -1) {

            bool concordant = isConcordant(newAlignment, alignments[pairTo]);
            if (bestNewPairScore >= bestPairScore) {
                concordantAlignmentExist = concordant;
            }
            if (bestNewPairScore > newAlignment->pairScore) {
                newAlignment->pairScore = bestNewPairScore;
                newAlignment->oppositePairAddresses.clear();
                newAlignment->oppositePairAddresses.push_back(alignments[pairTo]);
                newAlignment->pairToLocation = alignments[pairTo]->location;
                newAlignment->setConcordant(concordant);
            } else if (bestNewPairScore == newAlignment->pairScore) {
                newAlignment->oppositePairAddresses.push_back(alignments[pairTo]);
                newAlignment->pairToLocation = alignments[pairTo]->location;
                newAlignment->setConcordant(concordant);
            }
            if (bestNewPairScore > alignments[pairTo]->pairScore) {
                alignments[pairTo]->pairScore = bestNewPairScore;
                alignments[pairTo]->oppositePairAddresses.clear();
                alignments[pairTo]->oppositePairAddresses.push_back(newAlignment);
                alignments[pairTo]->pairToLocation = newAlignment->location;
                alignments[pairTo]->setConcordant(concordant);
            } else if (bestNewPairScore == alignments[pairTo]->pairScore) {
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
                uint64_t &multipleAlignment, uint64_t &repeatAlignment) {


        bool repeat = false;
        // find NH (nAlginment)
        int nAlignment = calculateNumAlignment(repeat);
        // update statistics
        if (nAlignment == 0) {
            unAligned++;
        } else {
            if (repeat && nAlignment>1 && !expandRepeat) {
                repeatAlignment++;
            } else {
                nonRepeatAlignment++;
                if (nAlignment == 1) {
                    uniqueAlignment++;
                } else {
                    multipleAlignment++;
                }
            }
        }
        // output
        bool primaryAlignment = true;
        int nOutput = 0;
        if (alignments.empty()) {
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

        assert ((nOutput == nAlignment) || (nAlignment == 0 && nOutput == 1));
        initialize();
    }

    void
    output_paired(BTString& o, uint64_t &unConcordant, uint64_t &uniqueConcordant, uint64_t &multipleConcordant, uint64_t &nonRepeatPairedAlignment,
                       uint64_t &uniqueDiscordant, uint64_t &unAlignedPairRead, uint64_t &alignedPairRead, uint64_t &uniqueAlignedPairRead, uint64_t &multipleAlignedPairRead,
                       uint64_t &concordantRepeat, uint64_t &disconcordantRepeatRead) {

        //int testNumBestPair = calculateNumBestPair();

        // update statistics
        if (!concordantAlignmentExist) {
            unConcordant++;
            int nAlignment[2] = {0};
            bool repeat[2] = {false};
            if (nBestPair == 0) {
                unAlignedPairRead += 2;
            } else if (nBestPair == 1) {
                uniqueDiscordant++;
            } else {
                for (int i = 0; i < alignments.size(); i++) {
                    if (alignments[i]->pairScore == bestPairScore || alignments[i]->repeat) {
                        nAlignment[alignments[i]->pairSegment]++;
                        if (alignments[i]->repeat) {
                            repeat[alignments[i]->pairSegment] = true;
                        }
                    }
                }
                for (int i = 0; i < 2; i++) {
                    if (repeat[i] && !expandRepeat) {
                        disconcordantRepeatRead++;
                    } else {
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
            }
        } else if (concordantAlignmentExist) {
            bool repeat = false;
            for (int i = 0; i < alignments.size(); i++) {
                if ((alignments[i]->pairScore == bestPairScore || alignments[i]->repeat) && alignments[i]->concordant == concordantAlignmentExist) {
                    if (alignments[i]->repeat && alignments[i]->nBestRepeat>1) {
                        repeat = true;
                    }
                }
            }
            if (repeat && !expandRepeat) {
                concordantRepeat++;
            } else {
                nonRepeatPairedAlignment++;
                if (nBestPair == 1) {
                    uniqueConcordant++;
                } else {
                    multipleConcordant++;
                }
            }
        }

        // output
        int nOutput = 0;
        bool primaryAlignment = true;
        if (alignments.empty() || nBestPair == 0) {
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


    void output(OutputQueue& oq_, size_t threadid_, uint64_t &unAligned, uint64_t &nonRepeatAlignment, uint64_t &uniqueAlignment,
                uint64_t &multipleAlignment, uint64_t &repeatAlignment,
                uint64_t &unConcordant, uint64_t &uniqueConcordant, uint64_t &multipleConcordant, uint64_t &nonRepeatPairedAlignment,
                uint64_t &uniqueDiscordant, uint64_t &unAlignedPairRead, uint64_t &alignedPairRead, uint64_t &uniqueAlignedPairRead, uint64_t &multipleAlignedPairRead,
                uint64_t &concordantRepeat, uint64_t &disconcordantRepeatRead) {

        while(working) {
            usleep(100);
        }
        BTString obuf_;
        OutputQueueMark qqm(oq_, obuf_, previousReadID, threadid_);
        if (paired) {
            output_paired(obuf_, unConcordant, uniqueConcordant, multipleConcordant, nonRepeatPairedAlignment,
                    uniqueDiscordant, unAlignedPairRead, alignedPairRead, uniqueAlignedPairRead, multipleAlignedPairRead,
                    concordantRepeat, disconcordantRepeatRead);
        } else {
            output_single(obuf_, unAligned, nonRepeatAlignment, uniqueAlignment, multipleAlignment, repeatAlignment);
        }
    }
};

#endif //HISAT2_TLA_H
