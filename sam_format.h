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

#ifndef HISAT2_SAM_FORMAT_H
#define HISAT2_SAM_FORMAT_H

#include <iostream>
using namespace std;


class SAMRecord {
public:
    SAMRecord() {};

    ~SAMRecord() {};


public:
    void parse_from_line(const string& linebuf);
    void update_tag(const string& name, const string& tagtype, const string& value);

    string m_qname; // Query name(read id)
    int m_flag;
    string m_rname; // Reference name
    int32_t m_pos;
    int m_mapq;
    string m_cigar;
    string m_rnext; // Reference name of the mate/next read
    int32_t m_pnext;   // Position of the mate/next read
    int32_t m_tlen; // observed Template length
    string m_seq;   // segment sequence
    string m_qual;  // Phred-scaled quality

    vector<string> m_options;   // optional fields

    friend ostream& operator<<(ostream& os, const SAMRecord& rec);
};
ostream& operator<<(ostream& os, const SAMRecord& rec);

#endif //HISAT2_SAM_FORMAT_H
