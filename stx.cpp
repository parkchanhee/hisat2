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

static inline pos_t calReadLengthCigar(const vector<CIGAR>& cigars)
{
    pos_t tot_len = 0;

    for (const auto& c: cigars) {
        if (c.op == 'M' || c.op == 'I' || c.op == 'S' || c.op == 'H') {
            tot_len += c.len;
        }
    }

    return tot_len;
}

void STXMap::loadFromMap(const string &fname)
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

            rec.stxChrname  = string(fields[0], 1);
            rec.chrname = fields[1];
            rec.gene = fields[2];
            rec.stxOffset = stoi(fields[3]);
            rec.stxLen = stoi(fields[4]);

            stxList.push_back(rec);
        } else {
            if (stxList.empty()) {
                cerr << "Invalid file format. Skip line: " << line << endl;
                continue;
            }

            STXRecord& current_rec = stxList.back();

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

    buildMaps();
}

void STXMap::loadFromTIDMap(const string &fname)
{
    ifstream fp(fname, ifstream::in);

    if (!fp.is_open()) {
        cerr << "Can't open file: " << fname << endl;
        return;
    }

    string line;

    STXRecord *pRec = NULL;
    while (getline(fp, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            continue;
        }

        if (line[0] == '>') {
            // header
            vector<string> fields;

            tokenize(line, delimiter, fields);
            assert(fields.size() >= 5);

            if (geneMap.find(fields[2]) != geneMap.end()) {
                pRec = geneMap[fields[2]];

                assert(pRec != NULL);

                pRec->bitWidth = stoi(fields[3]);
                pRec->numTIDs = stoi(fields[4]);
            } else {

                cerr << "No gene: " << fields[2] << endl;
                pRec = NULL;
            }

        } else {
            if (pRec == NULL) {
                cerr << "Invalid line. No headers: " << line << endl;
                continue;
            }
            vector<string> fields;

            tokenize(line, delimiter, fields);
            assert(fields.size() >= 2);

            BitVector bitVector(pRec->bitWidth);

            bitVector.load_from_hex(fields[0]);

            pRec->tidList.push_back(pair<BitVector, string>(bitVector, fields[1]));

#if 0
            cerr << "fields[0]: " << fields[0] << ", " << bitVector.to_hex() << ", tid: " << fields[1] << endl;

            if (fields[0] != bitVector.to_hex(true)) {
                cerr << "--------------Wrong------------------" << endl;
            }
#endif
        }
    }
}

void STXMap::saveToMap(const string &fname)
{
    ofstream fp(fname, ofstream::out);
    saveToMap(fp);
    fp.close();
}

void STXMap::saveToMap(ostream &os)
{
    for (const STXRecord& rec: stxList) {
        os  << ">" << rec.stxChrname << delimiter
            << rec.chrname << delimiter
            << rec.gene << delimiter
            << rec.stxOffset << delimiter
            << rec.stxLen << endl;

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

int STXMap::mapPosition(const string &stxChrname, pos_t stxPos, string &chrname, pos_t &posChr)
{
    pos_t offset = findOffset(stxChrname, stxPos);
    if (offset == INVALID_POS) {
        return -1;
    }

    // Get Record
    const record_map_t key = record_map_t(stxChrname, offset);
    const STXRecord& rec = *recordMap[key];

    pos_t new_pos;
    int ret = findPosition(rec, stxPos, new_pos);
    if (ret < 0) {
        cerr << "Can't find position: " << stxChrname << ":" << stxPos << endl;
        return 0;
    }

#ifdef DEBUGLOG
    cerr << stxChrname << ":" << stxPos << " -> " << rec.gene << ":" << offset << " -> " << rec.chrname << ":" << new_pos << endl;
#endif

    return 0;
}

int STXMap::mapPosition(const Position& inPosition, Position& outPosition)
{
    pos_t offset = findOffset(inPosition.chrname, inPosition.pos);
    if (offset == INVALID_POS) {
        return -1;
    }

    // Get Record
    const record_map_t key = record_map_t(inPosition.chrname, offset);
    const STXRecord& rec = *recordMap[key];

    return findPosition(rec, inPosition, outPosition, true, true);
}

int STXMap::mapPosition(SAMRecord& sam)
{
    //    Position old("22_tome", 2133882, "100M");
    Position old(sam.m_rname, sam.m_pos, sam.m_cigar);
    Position n;

    int ret = mapPosition(old, n);
    if (ret < 0) {
        return ret;
    }

    sam.m_rname = n.chrname;
    sam.m_pos = n.pos;
    sam.m_cigar = n.cigar;

    if (!n.gene.empty()) {
        // Append Gene
        sam.update_tag("GI", "Z", n.gene);
    }

    if (n.transcripts.size() > 0) {
        // Append Transcripts

        sam.update_tag("TI", "Z", n.transcripts[0]);

        if (n.transcripts.size() > 1) {
            stringstream ss;

            ss << n.transcripts[1];
            for (int i = 2; i < n.transcripts.size(); i++) {
                ss << "|" << n.transcripts[i];
            }

            sam.update_tag("TO", "Z", ss.str());
        }
    }
    return 0;
}

int STXMap::findPosition(const STXRecord &rec, pos_t pos, pos_t &posChr)
{
    // offset in SuperTranscript
    pos_t offset = pos - rec.stxOffset;
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

    posChr = new_pos;
    return 0;
}

int STXMap::findPosition(const STXRecord &rec, const Position &inPosition, Position &outPosition,
                          bool bCigar, bool bExonBitmap)
{
    pos_t offset_in_stx = inPosition.pos - rec.stxOffset;
    pos_t new_pos = INVALID_POS;

    int f_exon_idx = 0;

    if (offset_in_stx >= rec.stxLen) {
        cerr << "out of SuperTranscript" << endl;
        return -1;
    }

    for (; f_exon_idx < rec.exons.size(); f_exon_idx++) {
        const STXExon& exon = rec.exons[f_exon_idx];

        if (offset_in_stx < exon.len) {
            new_pos = exon.pos + offset_in_stx;
            break;
        } else {
            offset_in_stx -= exon.len;
        }
    }

    if (new_pos == INVALID_POS) {
        cerr << "Invalid Pos" << endl;
        return -2;
    }

    outPosition.chrname = rec.chrname;
    outPosition.pos = new_pos;

    if (!bCigar) {
        return 0;
    }

    if (inPosition.cigar.empty() || inPosition.cigar == "*") {
        return 0;
    }

    // Mapping CIGAR
//    cerr << "first exon idx: " << f_exon_idx << endl;
//    cerr << "offset_in_stx: " << offset_in_stx << endl;

    vector<CIGAR> cigars;
    cigarToList(inPosition.cigar, cigars);
    pos_t read_len = calReadLengthCigar(cigars);

    if (offset_in_stx + read_len > rec.stxLen) {
        fprintf(stderr, "Read length(%u) and offset(%u) is longer than stx length(%u)\n",
                read_len, offset_in_stx, rec.stxLen);
        return -1;
    }

    vector<int> exon_bit_list;
    exon_bit_list.push_back(f_exon_idx);

    vector<CIGAR> out_cigars;

    int e_idx = f_exon_idx;
    pos_t r_len = rec.exons[e_idx].len - (new_pos - rec.exons[e_idx].pos);

    for (const auto& c: cigars) {
        char c_op = c.op;
        pos_t c_len = c.len;

        if (c_op == 'S') {
            out_cigars.push_back(CIGAR(c_op, c_len));
            continue;
        }

        while (c_len > 0 && e_idx < rec.exons.size()) {
            if (c_len <= r_len) {
                r_len -= c_len;
                out_cigars.push_back(CIGAR(c_op, c_len));
                break;
            }

            c_len -= r_len;
            out_cigars.push_back(CIGAR(c_op, r_len));

            pos_t gap = INVALID_POS;
            if (e_idx == (rec.exons.size() - 1)) {
                // wrong alignment (across a transcript)
                gap = c_len;
                e_idx++;
            } else {
                gap = rec.exons[e_idx + 1].pos - (rec.exons[e_idx].pos + rec.exons[e_idx].len);
                e_idx++;
                r_len = rec.exons[e_idx].len;
            }

            if (c_op == 'M') {
                exon_bit_list.push_back(e_idx);
            }

            out_cigars.push_back(CIGAR('N', gap));
        }
    }

    if (exon_bit_list.back() == rec.exons.size()) {
        exon_bit_list.pop_back();
    }

    // merge cigars
    int ni = 0;
    for (int i = 1; i < out_cigars.size(); i++) {
        if (out_cigars[i].len == 0) {
            continue;
        } else if (out_cigars[i].op == 'S' && out_cigars[ni].op == 'N') {
            out_cigars[ni] = out_cigars[i];
        } else if (out_cigars[i].op == 'N' && out_cigars[ni].op == 'S') {
            continue;
        } else if (out_cigars[ni].op == out_cigars[i].op) {
            // merge
            out_cigars[ni].len += out_cigars[i].len;
        } else {
            ni++;
            out_cigars[ni] = out_cigars[i];
        }
    }
    out_cigars.resize(ni+1);

    outPosition.cigar = cigarListToString(out_cigars);
    outPosition.gene = rec.gene;


    if (!bExonBitmap) {
        return 0;
    }

    // Find TIDs
    //vector<string> tid_list;
    rec.findTIDs(exon_bit_list, outPosition.transcripts);

    return 0;
}
int STXMap::findPosition(const STXRecord &rec, pos_t pos, const vector<string> &cigars, pos_t &posChr,
                          vector<string> &cigarsChr)
{
    // offset in SuperTranscript
    pos_t offset = pos - rec.stxOffset;
    pos_t new_pos = INVALID_POS;

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

void STXMap::buildMaps()
{
    for (const auto& rec: stxList) {

        const string& stx_chrname = rec.stxChrname;

        offsetMap[stx_chrname].push_back(rec.stxOffset);

        record_map_t key = record_map_t(stx_chrname, rec.stxOffset);
        recordMap[key] = (STXRecord *)&rec;

        geneMap[rec.gene] = (STXRecord *)&rec;
    }

    for (auto& key_value: offsetMap) {
        sort(key_value.second.begin(), key_value.second.end());
    }

#ifdef DEBUGLOG
    cerr << "Total count of str_chrname: " << offsetMap.size() << endl;
    cerr << "Total number of record_map: " << recordMap.size() << endl;
#endif

}


pos_t STXMap::findOffset(const string &stxChrname, pos_t pos)
{
    const vector<pos_t>& offsets = offsetMap[stxChrname];
    if (offsets.size() == 0) {
        return INVALID_POS;
    }

    // offsets is sorted
    // find nearest offset of pos
    int idx = bisectRight(offsets, pos);

    if (idx < 0) {
        return INVALID_POS;
    }

    return offsets[idx];
}

int STXMap::bisectRight(const vector<pos_t> &list, pos_t pos)
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

int STXMap::cigarToList(const string &cigar, vector<CIGAR> &cigarList)
{
    const regex re("(\\d+)([MIDNSHPX=])");
    cigarList.resize(0);

    if (cigar.empty()) {
        return 0;
    }

    sregex_iterator itr(cigar.begin(), cigar.end(), re);
    sregex_iterator itr_end;

    for(; itr != itr_end; itr++) {
        std::smatch match = *itr;

        cigarList.push_back(CIGAR(match.str(2)[0], stoi(match.str(1))));
    }

    return cigarList.size();
}

string STXMap::cigarListToString(const vector<CIGAR> &cigarList)
{
    stringstream ss;

    for (const auto&c : cigarList) {
        ss << c.len << c.op;
    }

    return ss.str();
}

int STXRecord::findTIDs(const vector<int> &bit_list, vector<string> &out_tids) const
{
    BitVector bitVector(bitWidth);
    BitVector bitMask(bitWidth);

    bitVector.all_clear();
    bitMask.all_clear();

    // Mask
    for (int i = bit_list.front(); i <= bit_list.back(); i++) {
        bitMask[i] = 1;
    }

    for (const auto& i: bit_list) {
        bitVector[i] = 1;
    }

//    cerr << "bitMask: " << bitMask.to_str() << ", bitVector: " << bitVector.to_str() << endl;

    for (const auto& tidmap : tidList) {
        const BitVector& tid_bitvector = tidmap.first;

#if 0
        for (const auto& r : tidList) {
            cerr << r.first.to_str() << ", " << (tid_bitvector & bitMask).to_str() << ", " << r.second << endl;
        }
#endif

        if ((tid_bitvector & bitMask) == bitVector) {
            out_tids.push_back(tidmap.second);
        }
    }

    return 0;
}

ostream& operator<<(ostream& os, const STXExon& stx)
{
    os << stx.chrname << ":" << stx.pos << ":" << stx.len;
    return os;
}

ostream& operator<<(ostream& os, const STXRecord& rec)
{
    os  << rec.stxChrname << ":"
        << rec.chrname << ":"
        << rec.gene << ":"
        << rec.stxOffset << ":"
        << rec.stxLen;

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
    loadFromMap(fname);
    //loadFromTIDMap("genome.gt.tmap");
    loadFromTIDMap("22.gt.tmap");

    string tmp_str;
    pos_t tmp_pos;

    if (mapPosition("17_tome", 757220, tmp_str, tmp_pos) < 0) {
        cerr << "Can't find position" << endl;
    }

    if (mapPosition("2_tome", 358917, tmp_str, tmp_pos) < 0) {
        cerr << "Can't find position" << endl;
    }

    if (mapPosition("2_tome", 362166, tmp_str, tmp_pos) < 0 ) {
        cerr << "Can't find position" << endl;
    }

    if (mapPosition("22_tome", 0, tmp_str, tmp_pos) < 0) {
        cerr << "Can't find position" << endl;
    }
    if (mapPosition("22_tome", 3498419 + 5870 , tmp_str, tmp_pos) < 0) {
        cerr << "Can't find position" << endl;
    }
    if (mapPosition("1_tome", 4179 , tmp_str, tmp_pos) < 0) {
        cerr << "Can't find position" << endl;
    }

    const string cigar = "1M53N1187M355N12M";
    vector<CIGAR> cigars;
    cigarToList(cigar, cigars);

    for (const auto& c: cigars) {
        cerr << c.len << c.op << endl;
    }

    cerr << "Construct from list: " << cigarListToString(cigars) << endl;

    cigarToList("99M100", cigars);
    for (const auto& c: cigars) {
        cerr << c.len << c.op << endl;
    }

    pos_t tot_len = calReadLengthCigar(cigars);
    cerr << "read_len_cigar: " << tot_len << endl;
    cerr << "Construct from list: " <<  cigarListToString(cigars) << endl;


    //Position old("22_tome", 377316, "100M");
    Position old("22_tome", 2133882, "100M");
//    Position old("22_tome", 406553, "100M");
    Position n;

    int ret = mapPosition(old, n);
    cerr << "map_position: " << ret << endl;
    if (ret == 0) {
        cerr << "\t" << old.chrname << ":" << old.pos << ":" << old.cigar << endl;
        cerr << "\t" << n.chrname << ":" << n.pos << ":" << n.cigar << ":" << n.gene << ", " << "# tid: " << n.transcripts.size() << endl;

        for (const auto& t: n.transcripts) {
            cerr << "\t\t" << t << endl;
        }
    } else {
        cerr << "Can't find position" << endl;
    }


    return;
}
