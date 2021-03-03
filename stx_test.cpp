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

#include <iostream>
#include <strings.h>
#include <fstream>

#include "stx.h"
#include "bitvector.h"
#include "sam_format.h"

using namespace std;

string get_filename(const char *long_filename)
{
    if (long_filename == NULL) {
        return string();
    }

    const char *pos = rindex(long_filename, '/');
    if (pos == NULL) {
        return string(long_filename);
    }

    return string(pos + 1);
}

void show_usage(int argc, char *argv[])
{
    const string application_name = get_filename(argv[0]);
    cerr << application_name << " <Map Filename>" << endl;
    return;
}

void bitvector_test()
{
    BitVector bit(33);
    BitVector bit2(32);

    for (int i = 0; i < 12; i++) {
        bit[i] = 1;
        bit2[i] = 1;
    }

    if (bit == bit2) {
        cerr << "Same bitvector" << endl;
    } else {
        cerr << "Different" << endl;
    }

    float zz = bit[2];
    cout << "zz " << zz << endl;

    bit[0] = 0;
    bit[4] = 0;
    bit[7] = 0;

    if (bit == bit2) {
        cerr << "Same bitvector" << endl;
    } else {
        cerr << "Different" << endl;
    }


    bit[33] = 1;
    cerr << bit.to_str() << endl;
    cerr << bit.to_hex() << endl;


    BitVector bit_a(16);
    BitVector bit_b(16);
    BitVector bit_mask(16);

    bit_a.all_clear();
    bit_b.all_set();
    bit_mask.all_clear();


    bit_a[0] = 1; bit_a[1] = 0; bit_a[2] = 1;
    bit_mask[0] = 1;
    cerr << bit_a.to_hex() << " " << bit_b.to_hex() << " with mask: " << bit_mask.to_hex() << endl;

    if (bit_a.cmp_mask(bit_b, bit_mask)) {
        cerr << "cmp_mask true" << endl;
    }
}


void sam_parse(const string& fname)
{

    STXMap stxmap;
    stxmap.loadFromMap("22.gt.map");
    stxmap.loadFromTIDMap("22.gt.tmap");

    ifstream fp(fname, ifstream::in);

    if (!fp.is_open()) {
        cerr << "Can't open file: " << fname << endl;
        return;
    }

    string linebuf;

    while (getline(fp, linebuf)) {
        if (linebuf.empty()) {
            continue;
        }

        if (linebuf[0] == '@') {
            cerr << "Comment: " << linebuf << endl;
            continue;
        }

        SAMRecord sam;
        sam.parse_from_line(linebuf);

        cout << sam << endl;
        if (stxmap.mapPosition(sam) < 0) {
            cout << "Can't map" << endl;
        }
        cout << sam << endl;

    }

    fp.close();

}

int main(int argc, char *argv[])
{
    STXMap stxmap;

    cerr << "Hello World" << endl;

    if (argc < 2) {
        show_usage(argc, argv);
        return 0;
    }

//    stxmap.load_from_map(argv[1]);


//    string tmp_str;
//    pos_t tmp_pos;
//    stxmap.map_position("17_tome", 757220, tmp_str, tmp_pos);
//    stxmap.map_position("2_tome", 358917, tmp_str, tmp_pos);
//    stxmap.map_position("2_tome", 362166, tmp_str, tmp_pos);
//
//
//    stxmap.map_position("22_tome", 0, tmp_str, tmp_pos);
//    stxmap.map_position("22_tome", 3498419 + 5870 , tmp_str, tmp_pos);
//
//    stxmap.map_position("1_tome", 4179 , tmp_str, tmp_pos);

//    stxmap.save_to_map("aa.map");

//    stxmap.test(argv[1]);

//    bitvector_test();

    sam_parse(argv[1]);
    return 0;
}
