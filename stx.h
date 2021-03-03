/*
* Copyright 2021, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
*
* This file is part of HISAT 2.
*
* HISAT 2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* HISAT 2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __HISAT2_STX_H__
#define __HISAT2_STX_H__


#include <iostream>
//#include <algorithm>
#include <vector>
#include <map>
#include "bitvector.h"
#include "sam_format.h"

using namespace std;

typedef uint64_t pos_t;     // position type
const pos_t INVALID_POS = (pos_t)-1;

struct Position {
    Position() {};
    Position(const string& in_chrname, pos_t in_pos, const string& in_cigar)
            : chrname(in_chrname)
            , pos(in_pos)
            , cigar(in_cigar) {};

    string chrname;
    pos_t pos;
    string cigar;

    string gene;
    vector<string> transcripts;
};

struct STXExon {
public:
    string chrname;
    uint64_t pos;
    uint64_t len;
};
ostream& operator<<(ostream& os, const STXExon& stx);


struct STXRecord {
public:
    string stxChrname; //
    string chrname;
    string gene;
    uint32_t stxOffset; // offset in str_chr
    uint32_t stxLen;

    vector<STXExon> exons;

    int bitWidth;
    int numTIDs;
    vector<pair<BitVector, string>> tidList;

    int findTIDs(const vector<int>& bit_list, vector<string>& out_tids) const;

};
ostream& operator<<(ostream& os, const STXRecord& rec);

struct CIGAR {
    CIGAR() {};
    CIGAR(char in_op, pos_t in_len)
        : len(in_len), op(in_op) {};

    pos_t len;
    char op;
};

/**
 * SuperTranscript mapping information
 *
 * Can convert an alignment result to different coordination system
 *  - Genome
 *  - Gene
 *  - SuperTranscript
 *  - Transcript
 */
typedef pair<string, pos_t> record_map_t;

class STXMap {
public:
    STXMap();
    ~STXMap();

public:
    /**
     * Load mapping information from map file
     *
     * @param fname
     */
    void loadFromMap(const string& fname);
    void loadFromTIDMap(const string& fname);

    /**
     * Save mapping information to map file
     *
     * @param fname
     */
    void saveToMap(const string& fname);
    void saveToMap(ostream& os);

    /**
     * Map location to genomic location
     *
     * @param stx_chrname
     * @param stx_pos
     * @param chrname
     * @param pos_chr
     * @return
     */
    int mapPosition(const string& stxChrname, pos_t stxPos,
                     string& chrname, pos_t& posChr);

    int mapPosition(const Position& inPosition, Position& outPosition);
    int mapPosition(SAMRecord& sam);

private:

    // list of all SuperTranscript Records
    vector<STXRecord> stxList;

    // mapping pair(stx_chrname, stx_offset) to STXRecord
    map<record_map_t, STXRecord *> recordMap;
    // mapping gene name to STXRecord
    map<string, STXRecord *> geneMap;

    // list of stxOffsets per each stx_chrname
    map<string, vector<pos_t>> offsetMap;

private:

    /**
     * Make offsetMap, recordMap and geneMap from stxList
     */
    void buildMaps();

    /**
     * Find nearest offset of pos
     *
     * @param stx_chrname
     * @param pos
     * @return
     */
    pos_t findOffset(const string& stx_chrname, pos_t pos);

    /**
     * Convert Position to new
     *
     * @param rec
     * @param pos
     * @param pos_in_chr
     * @return
     */
    int findPosition(const STXRecord& rec, pos_t pos, pos_t& posChr);
    int findPosition(const STXRecord& rec, pos_t pos, const vector<string>& cigars, pos_t& posChr, vector<string>& cigarsChr);
    int findPosition(const STXRecord& rec, const Position& inPosition, Position& outPosition, bool bCigar = false, bool bExonBitmap = false);


public:
    void test(const string& fname);
    static int bisectRight(const vector<pos_t>& list, pos_t pos);
    static int cigarToList(const string& cigar, vector<CIGAR>& cigarList);
    static string cigarListToString(const vector<CIGAR>& cigarList);
};


#endif // __HISAT2_STX_H__