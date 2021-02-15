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

#include "sam_format.h"

using namespace std;

static const string delimiter = "\t";

ostream& operator<<(ostream& os, const SAMRecord& rec)
{
    string options_str;

    if (rec.m_options.size() > 0) {
        stringstream ss;
        for (const auto& t: rec.m_options) {
            ss << t << delimiter;
        }

        options_str = ss.str();
        if (options_str.size() > 0) {
            options_str.pop_back(); // remove last delimiter
        }

    }
    os  << rec.m_qname << delimiter
        << hex << rec.m_flag << dec << delimiter
        << rec.m_rname << delimiter
        << rec.m_pos << delimiter
        << rec.m_cigar
        ;
    if (options_str.size() > 0) {
        os << delimiter << options_str;
    }
    return os;
}

void SAMRecord::parse_from_line(const string &linebuf)
{

    vector<string> fields;

    tokenize(linebuf, delimiter, fields);
    assert(fields.size() > 10);


    m_qname = fields[0];
    m_flag = stoi(fields[1]);
    m_rname = fields[2];

    m_pos = stoi(fields[3]);
    m_mapq = stoi(fields[4]);

    m_cigar = fields[5];
    m_rnext = fields[6];

    m_pnext = stoi(fields[7]);
    m_tlen = stoi(fields[8]);

    m_seq = fields[9];
    m_qual = fields[10];

    m_options.assign(fields.begin() + 11, fields.end());
}


void SAMRecord::update_tag(const string &name, const string &tagtype, const string &value)
{
    bool found = false;
    string tag = name + ":" + tagtype + ":";

    for (int i = 0; i < m_options.size(); i++) {
        if (m_options[i].compare(0, tag.size(), tag) == 0) {
            m_options[i] = tag + value;
            found = true;
        }
    }

    if (!found) {
        // append
        m_options.push_back(tag + value);
    }

}