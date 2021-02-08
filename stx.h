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

using namespace std;


/**
 * SuperTranscript mapping information
 *
 * Can convert an alignment result to different coordination system
 *  - Genome
 *  - Gene
 *  - SuperTranscript
 *  - Transcript
 */
class STXMap {
public:
    STXMap();
    ~STXMap();
};


#endif // __HISAT2_STX_H__