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

#include <vector>
#include <string>
#include <iostream>

#include "bitvector.h"

using namespace std;


void BitVector::resize(size_t n_bit_len)
{
    data.resize(0);

    int data_len = (n_bit_len + bit_per_data - 1) / bit_per_data;
    data.resize(data_len);

    bit_len = n_bit_len;
}

int BitVector::get(size_t pos) const
{
    if (data.empty() || pos >= bit_len) {
        return -1;
    }

    size_t index = pos / bit_per_data;
    size_t offset = pos % bit_per_data;

    uint32_t t = data[index];

    return (t >> offset) & 0x1;
}

int BitVector::set(size_t pos, int value)
{
    if (data.empty() || pos >= bit_len) {
        return -1;
    }

    size_t index = pos / bit_per_data;
    size_t offset = pos % bit_per_data;

    uint32_t m = 1 << offset;

    if (value > 0) {
        // set
        data[index] |= m;
    } else {
        // clear
        data[index] &= ~m;
    }
    return value > 0;
}

bool BitVector::cmp_mask(const BitVector &b, const BitVector &m) const
{
    if (bit_len != b.bit_len || bit_len != m.bit_len) {
        return false;
    }

    for (int i = 0; i < bit_len; i++) {
        int v1 = get(i) & m.get(i);
        int v2 = b.get(i) & m.get(i);
        if (v1 != v2) {
            return false;
        }
    }

    return true;
}

void BitVector::all_set()
{
    for (int i = 0; i < bit_len; i++) {
        set(i, 1);
    }
}

void BitVector::all_clear()
{
    for (int i = 0; i < bit_len; i++) {
        set(i, 0);
    }
}

void BitVector::load_from_hex(const string &hex_str)
{
    //cerr << "hex_str:" << hex_str << endl;
    int bitpos = 0;

    string::const_reverse_iterator rit;
    for (rit = hex_str.rbegin(); rit != hex_str.rend(); rit++) {
        int val;
        char cv = *rit;
        if (isalpha(cv)) {
            val = toupper(*rit) - 'A' + 10;
        } else {
            val = *rit - '0';
        }

        if (val < 0 || val > 15) {
            cerr << "Invalid alnum: " << cv << endl;
            val = 0;
        }

        //cerr << *rit << ", " << val << endl;
        for (int i = 0; i < 4; i++) {
            if (bitpos >= bit_len) {
                return;
            }

            set(bitpos, (val >> i) & 0x1);

            bitpos++;
        }


    }
}


string BitVector::to_hex(bool bTrimLeadingZero) const
{
    char buf[32];
    string tstr;

    for (int i = 0; i < data.size() - 1; i++) {
        snprintf(buf, 32,"%08X", data[i]);
        tstr = string(buf) + tstr;
    }

    snprintf(buf, 32, "%08X", data.back());
    string msb_str(buf);

    int rlen = bit_len % bit_per_data;
    if (rlen) {
        rlen = (rlen + 3) / 4;

        tstr = msb_str.substr(msb_str.length() - rlen, rlen) + tstr;
    } else {
        tstr = msb_str + tstr;
    }

    if (bTrimLeadingZero) {
        int non_zero_index = 0;
        for (non_zero_index = 0; non_zero_index < tstr.size(); non_zero_index++) {
            if (tstr[non_zero_index] != '0') {
                break;
            }
        }

        if (non_zero_index > 0) {
            tstr = tstr.substr(non_zero_index, tstr.size() - non_zero_index);
        }

    }
    return tstr;
}

string BitVector::to_str() const
{
    string tstr;

    for(int i = bit_len - 1; i >= 0; i--) {
        tstr.push_back(get(i) ? '1':'0');
    }

    return tstr;
}

BitVector BitVector::operator&(const BitVector& b) const
{
    size_t new_bit_len = max(bit_len, b.bit_len);
    BitVector result(new_bit_len);

    result.all_clear();

    for (int i = 0; i < min(bit_len, b.bit_len); i++) {
        result.set(i, get(i) & b.get(i));
    }

    return result;
}
