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

using namespace std;

typedef uint64_t pos_t;     // position type

struct STXExon {
public:
    string chrname;
    uint64_t pos;
    uint64_t len;
};
ostream& operator<<(ostream& os, const STXExon& stx);


struct STXRecord {
public:
    string stx_chrname; //
    string chrname;
    string gene;
    uint32_t stx_offset; // offset in str_chr
    uint32_t stx_len;

    vector<STXExon> exons;
};
ostream& operator<<(ostream& os, const STXRecord& rec);

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
    void load_from_map(const string& fname);

    /**
     * Save mapping information to map file
     *
     * @param fname
     */
    void save_to_map(const string& fname);
    void save_to_map(ostream& os);

    /**
     * Map location to genomic location
     *
     * @param stx_chrname
     * @param stx_pos
     * @param chrname
     * @param pos_chr
     * @return
     */
    int map_position(const string& stx_chrname, const pos_t pos,
                     string& chrname, pos_t& pos_chr);

private:

    // list of all SuperTranscript Records
    vector<STXRecord> stx_list;


    // mapping pair(stx_chrname, stx_offset) to STXRecords
    map<record_map_t, STXRecord *> record_map;

    // list of stx_offsets per each stx_chrname
    map<string, vector<pos_t>> offset_map;

private:

    /**
     * Make offset_list from stx_list
     */
    void build_offset_list();

    /**
     * Make record_map from stx_list
     */
    void build_record_map();

    /**
     * Find nearest offset of pos
     *
     * @param stx_chrname
     * @param pos
     * @return
     */
    pos_t find_offset(const string& stx_chrname, pos_t pos);

    int find_position(const STXRecord& rec, pos_t pos, pos_t& pos_in_chr);
    int find_position(const STXRecord& rec, pos_t pos, const vector<string>& cigars, pos_t& pos_in_chr, vector<string>& cigars_in_chr);

public:
    void test(const string& fname);
    static int bisect_right(const vector<pos_t>& list, pos_t pos);
    static int cigar_to_list(const string& cigar, vector<string>& cigar_list);
};


#endif // __HISAT2_STX_H__