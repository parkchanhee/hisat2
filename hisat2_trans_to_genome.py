#!/usr/bin/env python3

import os, sys, math, random, re
from collections import defaultdict, Counter
from argparse import ArgumentParser, FileType
import pprint

cigar_re = re.compile('\d+\w')


def read_len_cigar(cigars):
    read_len = 0
    for cigar in cigars:
        length, cigar_op = cigar
        if cigar_op in "MISH":
            read_len += int(length)

    return read_len

def read_transinfo(info_fp):

    tbl = {}

    current_trans = None
    # parse
    for line in info_fp:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        if line.startswith('>'):
            line = line[1:].split('\t')

            current_tid = line[0]

            if not current_tid in tbl:
                #                   chr_name, strand, len, exons, gene_id
                tbl[current_tid] = [line[1], line[2], int(line[3]), list(), line[4]]
                current_trans = tbl[current_tid]

            else:
                print('Duplicated tid', current_tid)
                current_trans = tbl[current_tid]

        else:
            field = line.split('\t')

            for item in field:
                chr_name, genomic_position, exon_len = item.split(':')[0:3]
                current_trans[3].append([int(genomic_position), int(exon_len)])

    return tbl


def translate_position(trans_tbl, tid, tr_pos, cigars):

    trans = trans_tbl[tid]
    read_len = read_len_cigar(cigars)
   
    #print(tr_pos, read_len, trans[2], cigars)
    assert tr_pos + read_len <= trans[2]

    new_chr = trans[0]

    exons = trans[3]
    assert len(exons) >= 1


   
    e_idx = -1 
    for i, exon in enumerate(exons):
        # find first exon
        if tr_pos >= exon[1]:
            tr_pos -= exon[1]
        else:
            new_pos = exon[0] + tr_pos
            e_idx = i
            break

    assert e_idx >= 0

    new_cigar = list()

    r_len = exons[e_idx][1] - (new_pos - exons[e_idx][0])
    for cigar in cigars:
        c_len, c_op = cigar
        c_len = int(c_len)

        while c_len > 0 and e_idx < len(exons):
            if c_len <= r_len:
                r_len -= c_len
                new_cigar.append('{}{}'.format(c_len, c_op))
                break
            else:
                c_len -= r_len
                new_cigar.append('{}{}'.format(r_len, c_op))

                gap = exons[e_idx + 1][0] - (exons[e_idx][0] + exons[e_idx][1])
                new_cigar.append('{}{}'.format(gap, 'N'))

                e_idx +=1
                r_len = exons[e_idx][1]


    return new_chr, new_pos, ''.join(new_cigar)
    #
    #    new_chr, new_pos, new_cigar = translate_position(transinfo_table, tid, tr_pos, cigars)



def main(sam_file, transinfo_file):
    transinfo_table = read_transinfo(transinfo_file)

    #pprint.pprint(transinfo_table)

    for line in sam_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            print(line)
            continue

        fields = line.split('\t')

        #print(fields[0], fields[1])

        tid = fields[2]
        tr_pos = int(fields[3]) - 1
        cigar_str = fields[5]

        cigars = cigar_re.findall(cigar_str)
        cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]

        new_chr, new_pos, new_cigar = translate_position(transinfo_table, tid, tr_pos, cigars)
    
        print(fields[0], new_chr, new_pos + 1, new_cigar)
        

    return


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Translate the transcription-based postion to the genomic position')

    parser.add_argument('sam_file',
            nargs='?',
            type=FileType('r'),
            help='input SAM file')

    parser.add_argument('transinfo_file',
            nargs='?',
            type=FileType('r'),
            help='transcript position file')

    args = parser.parse_args()

    if not args.sam_file or not args.transinfo_file:
        parser.print_help()
        exit(1)

    main(args.sam_file, args.transinfo_file)

