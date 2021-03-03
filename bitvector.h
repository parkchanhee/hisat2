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


#ifndef __BITVECTOR_H__
#define __BITVECTOR_H__

#include <vector>
using namespace std;

const size_t bit_per_data = 32; // sizeof(uint32_t);

class BitVector {
public:
    BitVector() {
        resize(8);
    }

    BitVector(size_t N)
    {
        resize(N);
    };

    ~BitVector() {};


public:
    void resize(size_t n_bit_len);

    bool cmp_mask(const BitVector& b, const BitVector& m) const;

    void all_set();
    void all_clear();
    void load_from_hex(const string& hex_str);

private:
    vector<uint32_t> data;
    size_t bit_len;

    int get(size_t pos) const;
    int set(size_t pos, int value);

public:
    struct BitVector_Proxy {

        BitVector_Proxy(BitVector& inData, size_t inPos)
                : bitVector(inData)
                , pos(inPos)
        {};

        size_t pos;
        BitVector& bitVector;

        BitVector_Proxy& operator=(int value) {
            bitVector.set(pos, value);
            return *this;
        }
        operator int() const {
            return bitVector.get(pos);
        }
    };

public:
    BitVector_Proxy operator[](size_t pos)
    {
        return BitVector_Proxy(*this, pos);
    }

    bool operator==(const BitVector& b) const
    {
        if (this->bit_len != b.bit_len) {
            return false;
        }
        for (int i = bit_len - 1; i >= 0; i--) {
            if (this->get(i) != b.get(i)) {
                return false;
            }
        }
        return true;
    }
    bool operator!=(const BitVector& b) const
    {
        return !operator==(b);
    }

    BitVector operator&(const BitVector& b) const;

    string to_str() const;
    string to_hex(bool bTrimLeadingZero = false) const;

};



#endif //__BITVECTOR_H__
