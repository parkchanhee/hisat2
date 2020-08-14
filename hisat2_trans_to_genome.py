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
            chr_name = line[1]

            if len(line) < 3:
                strand = ''
                tran_gene_len = 100000000
                gene_id = ''
            else:
                strand = line[2]
                tran_gene_len = int(line[3])
                gene_id = line[4]


            if not current_tid in tbl:
                #                   chr_name, strand, len, exons, gene_id
                tbl[current_tid] = [chr_name, strand, tran_gene_len, list(), gene_id]
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
    new_pos = 0
    for i, exon in enumerate(exons):
        # find first exon
        if tr_pos >= exon[1]:
            tr_pos -= exon[1]
        else:
            new_pos = exon[0] + tr_pos
            e_idx = i
            break

    assert e_idx >= 0

    tmp_cigar = list()

    r_len = exons[e_idx][1] - (new_pos - exons[e_idx][0])
    for cigar in cigars:
        c_len, c_op = cigar
        c_len = int(c_len)

        while c_len > 0 and e_idx < len(exons):
            if c_len <= r_len:
                r_len -= c_len
                tmp_cigar.append([c_len, c_op])
                break
            else:
                c_len -= r_len
                tmp_cigar.append([r_len, c_op])

                gap = exons[e_idx + 1][0] - (exons[e_idx][0] + exons[e_idx][1])
                tmp_cigar.append([gap, 'N'])

                e_idx += 1
                r_len = exons[e_idx][1]


    # clean new_cigars
    ni = 0

    for i in range(1, len(tmp_cigar)):
        if tmp_cigar[i][0] == 0:
            pass
        elif tmp_cigar[i][1] == 'S' and tmp_cigar[ni][1] == 'N':
            # replace ni to i
            tmp_cigar[ni] = tmp_cigar[i]
        elif tmp_cigar[i][1] == 'N' and tmp_cigar[ni][1] == 'S':
            # goto next
            pass
        elif tmp_cigar[ni][1] == tmp_cigar[i][1]:
            # merge to ni
            tmp_cigar[ni][0] += tmp_cigar[i][0]
        else:
            ni += 1
            tmp_cigar[ni] = tmp_cigar[i]

    new_cigar = list()
    for i in range(ni+1):
        new_cigar.append('{}{}'.format(tmp_cigar[i][0], tmp_cigar[i][1]))

    return new_chr, new_pos, ''.join(new_cigar)
    #
    #    new_chr, new_pos, new_cigar = translate_position(transinfo_table, tid, tr_pos, cigars)


class OutputQueue:
    def __init__(self, fp=sys.stdout):
        self.fp = fp

        # read id
        self.current_rid = ''
        # key = [new_chr, new_pos, new_cigar], value = [count, samline_fields]
        self.alignments = {}

    def reset(self):
        self.current_rid = ''
        self.alignments = {}

        return

    """
    """
    def updateNH(self, fields, count):
        for i, tag in enumerate(fields):
            if isinstance(tag, str) and tag.startswith("NH:i:"):
                fields[i] = 'NH:i:{}'.format(str(count))
                break

        return

    def print_sam(self):

        count = len(self.alignments)

        for item, value in self.alignments.items():
            fields = value[1]

            flag = int(fields[1])
            mapq = 60
            if count == 1:
                mapq = 60
            else:
                mapq = 1

            # update NH
            self.updateNH(fields, count)

            assert self.current_rid == '' or self.current_rid == fields[0]

            print('\t'.join([fields[0], str(flag), item[0], str(item[1] + 1), str(mapq), item[2], *fields[6:]]), file=self.fp)


    """
    """
    def add(self, new_rid, new_chr, new_pos, new_cigar, sam_fields):
        new_key = (new_chr, new_pos, new_cigar)

        assert self.current_rid == '' or self.current_rid == new_rid

        if new_key in self.alignments:
            self.alignments[new_key][0] += 1
        else:
            self.alignments[new_key] = [1, sam_fields]

    """
    """
    def flush(self):
        if self.alignments:
            self.print_sam()

        self.reset()

    """
    """
    def push(self, new_rid, new_chr, new_pos, new_cigar, sam_fields):
        if self.current_rid != new_rid:
            self.flush()
            self.current_rid = new_rid

        self.add(new_rid, new_chr, new_pos, new_cigar, sam_fields)
        return

def main(sam_file, transinfo_file):
    transinfo_table = read_transinfo(transinfo_file)

    outq = OutputQueue()
    #pprint.pprint(transinfo_table)

    for line in sam_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            outq.flush()
            print(line)
            continue

        fields = line.split('\t')

        flag = int(fields[1])
        if flag & 0x04:
            # unmapped
            outq.flush()
            print(line)
            continue

        rid = fields[0]
        tid = fields[2]
        tr_pos = int(fields[3]) - 1
        cigar_str = fields[5]

        cigars = cigar_re.findall(cigar_str)
        cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]

        new_chr, new_pos, new_cigar = translate_position(transinfo_table, tid, tr_pos, cigars)

        outq.push(rid, new_chr, new_pos, new_cigar, fields)


        #print(fields[0], new_chr, new_pos + 1, new_cigar)
        
    outq.flush()

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

