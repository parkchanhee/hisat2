//
// Created by Yun (Leo) Zhang on 1/10/20.
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

    RepeatPosition(long long int &inputLocation, string &inputChromosome):location(inputLocation), chromosome(inputChromosome) {
        pairScore = numeric_limits<int>::min();
        pairToLocation = -1;
        concordant = false;
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

    MappingPosition (long long int inputLocation, string inputChromosome, int inputAS_score, int inputReadSegment=0,
                     int inputPairScore=numeric_limits<int>::min(), long long int inputPairToLocation = -1){
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


    int calculatePairScore(int &AS1, int &AS2, string &chromosome1, string &chromosome2, long long int &location1, long long int &location2) {
        if (chromosome1 == chromosome2) {
            return (10*AS1+10*AS2-abs(location1-location2));
        } else {
            return (numeric_limits<int>::min()/2);
        }
    }

    int getPairScore(int &inputAS, string &inputChromosome, long long int &inputLocation) {
        int score = calculatePairScore(AS, inputAS, chromosome, inputChromosome, location, inputLocation);
        if (score > pairScore) {
            pairScore = score;
            pairToLocation = inputLocation;
        }
        return score;
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

    void appendPair (long long int location, string chromosome, int AS, int pairSegment=0) {

        int score;
        int index;
        int bestScoreThisPair = numeric_limits<int>::min();
        long long int bestPairToLocation;

        for (int i = 0; i < positions.size(); i++) {
            if (pairSegment == positions[i].pairSegment) {
                continue;
            } else {
                score = positions[i].getPairScore(AS, chromosome, location);
                if (score > bestScoreThisPair) {
                    bestScoreThisPair = score;
                    bestPairToLocation = positions[i].location;
                }
            }
        }

        if (bestScoreThisPair > bestPairScore) {
            bestPairScore = bestScoreThisPair;
            nPair = 1;
        } else if (bestScoreThisPair == bestPairScore) {
            nPair++;
        }

        positions.push_back(MappingPosition(location, chromosome, AS, pairSegment, bestScoreThisPair, bestPairToLocation));
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
    long long int* refLocations = nullptr;
    char* locationPointer;
    bool planA;
    bool paired;
    bool forward;
    bool primaryAlignment;
    bool mapped;
    bool concordant;
    int pairSegment;
    struct ht2_repeat_expand_result *repeatResult = nullptr;
    int pairScore; // to identify the better pair

    // for repeatAlignment only
    bool repeat;
    bool pairToRepeat;
    MappingPositions repeatPositions; // only have chromosme and location informations
    int bestAS;
    int nBestRepeat;
    Alignment* oppositePairAddress;

    void initialize() {
        readName.clear();
        flag = -1;
        chromosomeName.clear();
        location = -1;
        MAPQ.clear();
        cigarString.clear();
        pairToChromosome.clear();
        pairToLocation = -1;
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
        oppositePairAddress = nullptr;
        repeatPositions.initialize();
        if (refLocations != nullptr) {
            free(refLocations);
            refLocations = nullptr;
        }
        if (repeatResult != nullptr) {
            free(repeatResult);
            repeatResult = nullptr;
        }
    }

    Alignment() {
        initialize();
    }

    ~Alignment() {
        if (refLocations != nullptr) {
            free(refLocations);
        }
        if (repeatResult != nullptr) {
            free(repeatResult);
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
    }

    int RA_Map(char &base) {
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
        // update extra RA tag
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

        ht2_error_t err = ht2_repeat_expand((planA == true)?repeatHandles[0]:repeatHandles[1],
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
            if (!constructMD(readSequence_string, newMD, nIncorrectMatch)) {
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

    bool constructMD(string &readSequence_string, string &newMD_String, int &nIncorrectMatch) {
        // construct new MD string, this function if for repeat read
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
        if (concordant) {
            o.append("=");
            o.append('\t');
        } else {
            o.append(pairToChromosome.toZBuf());
            o.append('\t');
        }

        // pair to location
        itoa10<int>(pairToLocation, buf);
        o.append(buf);
        o.append('\t');

        // pairing distance
        if (paired) {
            itoa10<int>(pairToLocation - location, buf);
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
        outputFlags(o);

        o.append('\n');
    }

    void outputRepeatOppositePair(BTString& o, Alignment* oppositeAlignment, string &chromosome, long long int &location, bool primary) {
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
                    o.append(oppositeAlignment->pairToChromosome.toZBuf());
                    o.append('\t');

                    // pair to location
                    itoa10<int>(currentPosition->repeatPositions[j].pairToLocation, buf);
                    o.append(buf);
                    o.append('\t');

                    // pairing distance
                    itoa10<int>(currentPosition->repeatPositions[j].pairToLocation-currentPosition->repeatPositions[j].location, buf);
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
                if (!paired && (currentPosition->AS == bestAS) || (paired && (currentPosition->repeatPositions[j].pairScore == pairScore))) {
                    nOutput++;

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
                    o.append(pairToChromosome.toZBuf());
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

                    if (oppositePairAddress != nullptr && oppositePairAddress->repeat && paired) {
                        outputRepeatOppositePair(o, oppositePairAddress, currentPosition->repeatPositions[j].chromosome, currentPosition->repeatPositions[j].pairToLocation, (nOutput==1 && primaryAlignment));
                    }
                }
            }
        }
        return nOutput;
    }

    int output(BTString& o) {
        int nOutput;
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
    vector<Alignment*> alignments;
    queue<Alignment*> freeAlignments;

    TReadId previousReadID;
    int bestAS; // for single-end read selection
    int bestPairScore; // for paired-end read selection;
    int nBestPair;
    bool concordantAlignmentExist;
    MappingPositions existPositions;

    BTString readName[2];
    BTString readSequence[2];
    BTString qualityScore[2];

    bool working;
    bool paired;

    void initialize() {
        //alignments.clear();
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

    int calculatePairScore(Alignment *alignment1, Alignment *alignment2, int &nPair) {
        int pairScore = numeric_limits<int>::min();

        nPair = 0;
        if (alignment1->pairSegment == alignment2->pairSegment){
            pairScore = numeric_limits<int>::min();
        } else if (!alignment1->mapped || !alignment2->mapped) {
            pairScore = numeric_limits<int>::min()/2;
            if (!alignment1->mapped && !alignment2->mapped) {
                nPair = 0;
                pairScore -= 1;
            } else {
                if (alignment1->chromosomeName == alignment2->chromosomeName) {
                    pairScore += 1;
                }
                nPair = 1;
            }

            if (pairScore > alignment1->pairScore) {
                alignment1->pairScore = pairScore;
            }
            if (pairScore > alignment2->pairScore) {
                alignment2->pairScore = pairScore;
            }

        } else if ((alignment1->chromosomeName != alignment2->chromosomeName && ((alignment1->chromosomeName != alignment2->pairToChromosome) && (alignment1->pairToChromosome != alignment2->chromosomeName))) ||
                    ((alignment1->location != alignment2->pairToLocation) && (alignment1->pairToLocation != alignment2->location))){
            //pairScore = numeric_limits<int>::min()/2 + 2;
            if (pairScore > alignment1->pairScore) {
                alignment1->pairScore = pairScore;
            }
            if (pairScore > alignment2->pairScore) {
                alignment2->pairScore = pairScore;
            }
            //nPair = 1;
        } else if ((!alignment1->repeat && !alignment2->repeat) || !expandRepeat){
            pairScore = calculatePairScore(alignment1->location, alignment1->AS, alignment2->location, alignment2->AS);
            pairScore += alignment1->concordant * 1000;
            if (pairScore > alignment1->pairScore) {
                alignment1->pairScore = pairScore;
                alignment1->pairToLocation = alignment2->location;
            }
            if (pairScore > alignment2->pairScore) {
                alignment2->pairScore = pairScore;
                alignment2->pairToLocation = alignment1->location;
            }
            nPair = 1;
        } else if (alignment1->repeat && !alignment2->repeat) {
            int tmpPairScore;

            int AS2 = alignment2->AS;
            long long int location2 = alignment2->location;
            string chrmosome2 = alignment2->chromosomeName.toZBuf();

            for (int i = 0; i < alignment1->repeatPositions.positions.size(); i++) {
                MappingPosition *repeatPositions = &alignment1->repeatPositions.positions[i];
                for (int j = 0; j < repeatPositions->repeatPositions.size(); j++) {
                    RepeatPosition *repeatPosition = &repeatPositions->repeatPositions[j];
                    if(repeatPosition->chromosome != chrmosome2) {
                        continue;
                    }
                    tmpPairScore = calculatePairScore(repeatPosition->location, repeatPositions->AS, location2, AS2);
                    tmpPairScore += alignment1->concordant * 1000;
                    if (tmpPairScore >= alignment1->pairScore) {
                        alignment1->pairScore = tmpPairScore;
                        repeatPositions->pairScore = tmpPairScore;
                        repeatPosition->pairScore = tmpPairScore;
                        repeatPosition->pairToLocation = location2;
                    }
                    if (tmpPairScore >= alignment2->pairScore) {
                        alignment2->pairScore = tmpPairScore;
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
            int tmpPairScore;

            int AS1 = alignment1->AS;
            long long int location1 = alignment1->location;
            string chromosome1 = alignment1->chromosomeName.toZBuf();

            for (int i = 0; i < alignment2->repeatPositions.positions.size(); i++) {
                MappingPosition *repeatPositions = &alignment2->repeatPositions.positions[i];
                for (int j = 0; j < repeatPositions->repeatPositions.size(); j++) {
                    RepeatPosition *repeatPosition = &repeatPositions->repeatPositions[j];
                    if(chromosome1 != repeatPosition->chromosome) {
                        continue;
                    }
                    tmpPairScore = calculatePairScore(repeatPosition->location, repeatPositions->AS, location1, AS1);
                    tmpPairScore += alignment1->concordant * 1000;
                    if (tmpPairScore >= alignment1->pairScore) {
                        alignment1->pairScore = tmpPairScore;
                        alignment1->pairToLocation = repeatPosition->location;
                    }
                    if (tmpPairScore >= alignment2->pairScore) {
                        alignment2->pairScore = tmpPairScore;
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
                            tmpPairScore = calculatePairScore(repeatPosition1->location, repeatPositions1->AS, repeatPosition2->location, repeatPositions2->AS);
                            tmpPairScore += alignment1->concordant * 1000;
                            if (tmpPairScore >= alignment1->pairScore) {
                                alignment1->pairScore = tmpPairScore;
                                repeatPositions1->pairScore = tmpPairScore;
                                repeatPosition1->pairScore = tmpPairScore;
                                repeatPosition1->pairToLocation = repeatPosition2->location;
                            }
                            if (tmpPairScore >= alignment2->pairScore) {
                                alignment2->pairScore = tmpPairScore;
                                repeatPositions2->pairScore = tmpPairScore;
                                repeatPosition2->pairScore = tmpPairScore;
                                repeatPosition2->pairToLocation = repeatPosition1->location;
                            }
                            if (tmpPairScore > pairScore) {
                                pairScore = tmpPairScore;
                                nPair = 1;
                            } else if (tmpPairScore == pairScore) {
                                nPair++;
                            }
                        }
                    }
                }
            }
        }
        return pairScore;
    }


    int calculatePairScore(long long int location1, int AS1, long long int location2, int AS2) {
        int score = 100*AS1 + 100*AS2;
        int distance = abs(location1 - location2);
        if (distance >= 10000) {
            score -= distance;
        }
        if (distance > 500000) {
            return numeric_limits<int>::min();
        }
        return score;
    }

    // add a extra tag if the read is aligned on repeat index. this tag will show which repeat index is used.
    // XR:Z:TC means the read aligned on TC conversion repeat index.
    void addRepeatIndexTag(Alignment *newAlignment) {
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

        working = true;

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
            for (int i = 0; i < alignments.size(); i++) {
                alignments[i]->initialize();
                freeAlignments.push(alignments[i]);
            }
            alignments.clear();
            bestAS = newAlignment->AS;
            alignments.push_back(newAlignment);
        } else if (newAlignmentAS == bestAS) {
            if (alignmentExist(newAlignment)) {
                newAlignment->initialize();
                freeAlignments.push(newAlignment);
            } else {
                alignments.push_back(newAlignment);
            }
        } else {
            newAlignment->initialize();
            freeAlignments.push(newAlignment);
        }
        working = false;
    }

    void addNewAlignment_paired(Alignment *newAlignment) {
        working = true;

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
                // find it's pair
                if ((alignments[i]->location != newAlignment->pairToLocation) || (alignments[i]->chromosomeName != newAlignment->chromosomeName)) {
                    continue;
                }
                if (alignments[i]->pairSegment == newAlignment->pairSegment) {
                    continue;
                }
                int nPair = 0;
                int newPairScore = calculatePairScore(newAlignment, alignments[i], nPair);
                if (newPairScore > bestNewPairScore) {
                    bestNewPairScore = newPairScore;
                    tmp_nBestPair = nPair;
                    pairTo = i;
                    if (newAlignment->concordant) {
                        concordantAlignmentExist = true;
                    }
                }
            }
        }
        if (bestNewPairScore > bestPairScore) {
            bestPairScore = bestNewPairScore;
            nBestPair = tmp_nBestPair;
        } else if (bestNewPairScore == bestPairScore) {
            nBestPair += tmp_nBestPair;
        }
        if (pairTo > -1) {
            if (bestNewPairScore > newAlignment->pairScore) {
                newAlignment->oppositePairAddress = alignments[pairTo];
                newAlignment->pairToLocation = alignments[pairTo]->location;
            }
            if (bestNewPairScore > alignments[pairTo]->pairScore) {
                alignments[pairTo]->oppositePairAddress = newAlignment;
                alignments[pairTo]->pairToLocation = newAlignment->location;
            }

            if (newAlignment->pairSegment == 0) {
                //alignments.insert(alignments.begin()+pairTo, newAlignment);
                alignments.push_back(newAlignment);
            } else {
                alignments.insert(alignments.begin()+pairTo+1, newAlignment);
            }
        } else {
            alignments.push_back(newAlignment);
        }

        working = false;
    }

    void outputUnAlignmentRead(BTString& o) {
        if (paired) {
            for (int i = 0; i < 2; i++) {
                o.append(readName[i].toZBuf());
                o.append("\t4\t*\t0\t0\t*\t*\t0\t0\t");
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

        int nAlignment = 0;
        bool repeat = false;
        // find NH (nAlginment)
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

    void output_paired(BTString& o, uint64_t &unConcordant, uint64_t &uniqueConcordant, uint64_t &multipleConcordant, uint64_t &nonRepeatPairedAlignment,
                       uint64_t &uniqueDiscordant, uint64_t &unAlignedPairRead, uint64_t &alignedPairRead, uint64_t &uniqueAlignedPairRead, uint64_t &multipleAlignedPairRead,
                       uint64_t &concordantRepeat, uint64_t &disconcordantRepeatRead) {
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
                    alignments[i]->output(o);
                    if (alignments[i]->pairSegment == 1) {
                        primaryAlignment = false;
                    }
                }
            }
        }
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
