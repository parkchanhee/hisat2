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
    void resize(size_t n_bit_len)
    {
        data.resize(0);

        int data_len = (n_bit_len + bit_per_data - 1) / bit_per_data;
        data.resize(data_len);

        bit_len = n_bit_len;
    }

    int get(size_t pos) const
    {
        if (data.empty() || pos >= bit_len) {
            return -1;
        }

        size_t index = pos / bit_per_data;
        size_t offset = pos % bit_per_data;

        uint32_t t = data[index];

        return (t >> offset) & 0x1;
    }

    int set(size_t pos, int value)
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

    bool cmp_mask(const BitVector& b, const BitVector& m) const
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

    void all_set() {
        for (int i = 0; i < bit_len; i++) {
            set(i, 1);
        }
    }
    void all_clear() {
        for (int i = 0; i < bit_len; i++) {
            set(i, 0);
        }
    }

private:
    vector<uint32_t> data;
    size_t bit_len;

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

    string to_str()
    {
        string tstr;

        for(int i = bit_len - 1; i >= 0; i--) {
            tstr.push_back(get(i) ? '1':'0');
        }

        return tstr;
    }

    string to_hex()
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
        }
        return tstr;
    }

};



#endif //__BITVECTOR_H__
