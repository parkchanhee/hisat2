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

#include "stx.h"

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

    stxmap.test(argv[1]);

    return 0;
}
