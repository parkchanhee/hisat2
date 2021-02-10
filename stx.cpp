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

#include <cassert>

#include <iostream>
#include <fstream>

#include <algorithm>
#include <regex>

#include "tokenize.h"
#include "stx.h"

using namespace std;

static const string delimiter = "\t";
static const string delimiter_exon = ":";

STXMap::STXMap()
{
}

STXMap::~STXMap()
{
}

void STXMap::load_from_map(const string &fname)
{
    ifstream fp(fname, ifstream::in);

    if (!fp.is_open()) {
        cerr << "Can't open file: " << fname << endl;
        return;
    }

    string line;

    while (getline(fp, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            // Header Line
            // stx_chrname chrname gene offset_in_stxchr stx_len
            vector<string> fields;
            STXRecord rec;

            tokenize(line, delimiter, fields);
            assert(fields.size() >= 5);

            rec.stx_chrname  = string(fields[0], 1);
            rec.chrname = fields[1];
            rec.gene = fields[2];
            rec.stx_offset = stoi(fields[3]);
            rec.stx_len = stoi(fields[4]);

            stx_list.push_back(rec);
        } else {
            if (stx_list.empty()) {
                cerr << "Invalid file format. Skip line: " << line << endl;
                continue;
            }

            STXRecord& current_rec = stx_list.back();

            vector<string> exon_maps;
            tokenize(line, delimiter, exon_maps);
            for (const string& ex: exon_maps) {
                vector<string> fields;

                tokenize(ex, delimiter_exon, fields);

                if (fields.size() != 3) {
                    cerr << "Invalid line. Skip item: " << ex << endl;
                    continue;
                }

                STXExon stx_exon;
                stx_exon.chrname = fields[0];
                stx_exon.pos = stol(fields[1]);
                stx_exon.len = stol(fields[2]);

                assert(stx_exon.chrname == current_rec.chrname);

                current_rec.exons.push_back(stx_exon);
            }
        }
    }

    fp.close();

    build_offset_list();
    build_record_map();
}

void STXMap::save_to_map(const string &fname)
{
    ofstream fp(fname, ofstream::out);
    save_to_map(fp);
    fp.close();
}

void STXMap::save_to_map(ostream &os)
{
    for (const STXRecord& rec: stx_list) {
        os  << ">" << rec.stx_chrname << delimiter
            << rec.chrname << delimiter
            << rec.gene << delimiter
            << rec.stx_offset << delimiter
            << rec.stx_len << endl;

        int count = 0;
        for (const STXExon& exon: rec.exons) {
            if (count > 0) {
                if (count % 4) {
                    os << delimiter;
                } else {
                    os << endl;
                }
            }

            os  << exon.chrname << delimiter_exon
                << exon.pos << delimiter_exon
                << exon.len;

            count++;
        }

        if (count > 0) {
            os << endl;
        }
    }
}

int STXMap::map_position(const string &stx_chrname, const pos_t pos, string &chrname, pos_t &pos_chr)
{
    pos_t offset = find_offset(stx_chrname, pos);

    // Get Record
    const record_map_t key = record_map_t(stx_chrname, offset);
    const STXRecord& rec = *record_map[key];

    pos_t new_pos;
    int ret = find_position(rec, pos, new_pos);
    if (ret < 0) {
        cerr << "Can't find position: " << stx_chrname << ":" << pos << endl;
        return 0;
    }

    cerr << stx_chrname << ":" << pos << " -> " << rec.gene << ":" << offset << " -> " << rec.chrname << ":" << new_pos << endl;

    return 0;
}

int STXMap::find_position(const STXRecord &rec, pos_t pos, pos_t &pos_in_chr)
{
    // offset in SuperTranscript
    pos_t offset = pos - rec.stx_offset;
    pos_t new_pos = 0;

    for(int i = 0; i < rec.exons.size(); i++) {
        const STXExon& exon = rec.exons[i];

        if (offset < exon.len) {
            new_pos = exon.pos + offset;
            break;
        } else {
            offset -= exon.len;
        }
    }

    pos_in_chr = new_pos;
    return 0;
}

int STXMap::find_position(const STXRecord &rec, pos_t pos, const vector<string> &cigars, pos_t &pos_in_chr,
                          vector<string> &cigars_in_chr)
{
    // offset in SuperTranscript
    pos_t offset = pos - rec.stx_offset;
    pos_t new_pos = (pos_t)-1;

    for(int i = 0; i < rec.exons.size(); i++) {
        const STXExon& exon = rec.exons[i];

        if (offset < exon.len) {
            new_pos = exon.pos + offset;
            break;
        } else {
            offset -= exon.len;
        }
    }

    if (new_pos == (pos_t)-1) {
        cerr << "Can't find new_pos" << endl;
        return -1;
    }

    return 0;
}

void STXMap::build_offset_list()
{
    for (const auto& i: stx_list) {
        const string& stx_chrname = i.stx_chrname;
        offset_map[stx_chrname].push_back(i.stx_offset);
    }

    for (auto& key_value: offset_map) {
        sort(key_value.second.begin(), key_value.second.end());
    }

    cerr << "Total count of str_chrname: " << offset_map.size() << endl;

}

void STXMap::build_record_map()
{
    for (const auto& rec: stx_list) {
        record_map_t key = record_map_t(rec.stx_chrname, rec.stx_offset);
        record_map[key] = (STXRecord *)&rec;
    }

    cerr << "Total number of record_map: " << record_map.size() << endl;
}

pos_t STXMap::find_offset(const string &stx_chrname, pos_t pos)
{
    const vector<pos_t>& offsets = offset_map[stx_chrname];
    if (offsets.size() == 0) {
        return 0;
    }

    // offsets is sorted
    // find nearest offset of pos
    int idx = bisect_right(offsets, pos);

    if (idx < 0) {
        return 0;
    }

    return offsets[idx];
}

int STXMap::bisect_right(const vector<pos_t> &list, pos_t pos)
{
    int lo = 0;
    int hi = list.size();


    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (pos < list[mid]) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }

    assert(lo > 0);
    return lo - 1;

}

int STXMap::cigar_to_list(const string &cigar, vector<string> &cigar_list)
{
    const regex re("\\d+\\w");
    cigar_list.resize(0);

    if (cigar.empty()) {
        return 0;
    }

    sregex_token_iterator itr(cigar.begin(), cigar.end(), re);
    sregex_token_iterator itr_end;

    for(; itr != itr_end; itr++) {
        cigar_list.push_back(*itr);
    }

    return cigar_list.size();
}

ostream& operator<<(ostream& os, const STXExon& stx)
{
    os << stx.chrname << ":" << stx.pos << ":" << stx.len;
    return os;
}

ostream& operator<<(ostream& os, const STXRecord& rec)
{
    os  << rec.stx_chrname << ":"
        << rec.chrname << ":"
        << rec.gene << ":"
        << rec.stx_offset << ":"
        << rec.stx_len;

    if (rec.exons.size() > 0) {
        os << endl;

        int count = 0;
        for (const auto& i : rec.exons) {
            os << i;
            count++;
            if (count %4 == 0) {
                os << endl;
            } else {
                os << " ";
            }
        }
        if (count > 0 && count % 4) {
            os << endl;
        }
    }

    return os;
}


void STXMap::test(const string& fname)
{
    load_from_map(fname);

    string tmp_str;
    pos_t tmp_pos;

    map_position("17_tome", 757220, tmp_str, tmp_pos);
    map_position("2_tome", 358917, tmp_str, tmp_pos);
    map_position("2_tome", 362166, tmp_str, tmp_pos);


    map_position("22_tome", 0, tmp_str, tmp_pos);
    map_position("22_tome", 3498419 + 5870 , tmp_str, tmp_pos);
    map_position("1_tome", 4179 , tmp_str, tmp_pos);


    vector<string> cigars;
    cigar_to_list("1M53N87M355N12M", cigars);

    for (const auto& i: cigars) {
        cerr << i << endl;
    }




    return;
}