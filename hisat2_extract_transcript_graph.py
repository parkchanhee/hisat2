#!/usr/bin/env python3

#
# Copyright 2020, Chanhee Park <parkchanhee@gmail.com>, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#

from sys import stderr, exit
from collections import defaultdict, Counter
from argparse import ArgumentParser, FileType
import pprint
import bisect

bDebug = False
bVerbose = False

"""
"""
def read_genome(genome_file, chr_filter = None):
    chr_dic = {}
    
    chr_name, sequence = "", ""
    for line in genome_file:
        if line[0] == ">":
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line[:-1]

    if chr_name and sequence:
        chr_dic[chr_name] = sequence

    #chr_filter = [str(x) for x in list(range(1, 23)) + ['X', 'Y']]
    #chr_filter = None

    if chr_filter:
        for chr_id, chr_seq in chr_dic.items():
            if not chr_id in chr_filter: 
                chr_dic.pop(chr_id, None)
    
    return chr_dic


"""
"""
def write_fasta(fp, ref_id, seq, width=60):
    print('>{}'.format(ref_id), file=fp)

    for i in range(0, len(seq), width):
        print(seq[i:i+width], file=fp)
    
    return


"""
"""
def load_transcript(genome_seq, gtf_file):
    genes = defaultdict(list)
    transcripts = {}

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
        if bVerbose:
            line_org = str(line)

        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()
        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            if bVerbose:
                print("Warning: Can't parse line:", line_org)

            continue
        if chrom not in genome_seq:
            if bVerbse:
                print('chr {} is not in genome'.format(chrom))

            continue
        
        # Zero-based offset
        left, right = int(left) - 1, int(right) - 1
        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';')[:-1]:
            attr, _, val = attr.strip().partition(' ')
            values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            print('wrong line', line)
            continue

        transcript_id = values_dict['transcript_id']
        gene_id = values_dict['gene_id']

        if transcript_id not in transcripts:
            transcripts[transcript_id] = [chrom, strand, [[left, right]], gene_id]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            transcripts[transcript_id][2].append([left, right])

    return genes, transcripts



"""
"""
def read_transcript(genome_seq, gtf_file, min_transcript_len = 0):
    genes, transcripts = load_transcript(genome_seq, gtf_file)

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chr, strand, exons, gene_id] in transcripts.items():
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tmp_exons[-1][1] = exons[i][1]
                else:
                    tmp_exons.append(exons[i])
            transcripts[tran] = [chr, strand, tmp_exons, gene_id]

    tmp_transcripts = {}
    for tran, [chr, strand, exons, gene_id] in transcripts.items():
        exon_lens = [e[1] - e[0] + 1 for e in exons]
        transcript_len = sum(exon_lens)
        if transcript_len >= min_transcript_len:
            tmp_transcripts[tran] = [chr, strand, transcript_len, exons, gene_id]

    transcripts = tmp_transcripts

    return genes, transcripts




"""
"""
def merge_range(A, B):
    assert A[0] <= B[0]

    if (A[1]+1) < B[0]:
        return None

    return [A[0], max(A[1], B[1])]

"""
"""
def make_consensus_exons(exons_list):

    assert len(exons_list) >= 1

    exons_list = sorted(exons_list)

    assert exons_list[0][0] <= exons_list[0][1]
    cons_exon_list = [exons_list[0]]

    for i in range(1, len(exons_list)):
        last_exon = cons_exon_list[-1]
        curr_exon = exons_list[i]

        assert curr_exon[0] <= curr_exon[1]

        # check two ranges
        m = merge_range(last_exon, curr_exon)
        if m is None:
            cons_exon_list.append([*curr_exon])
        else:
            cons_exon_list[-1] = m


    accum_len = 0
    for e in cons_exon_list:
        e.append(accum_len)
        accum_len += e[1] - e[0] + 1

    return cons_exon_list


"""
"""
def write_transcripts_seq(fp, gene_id, chr_name, genome_seq, exon_list):
    chr_seq = genome_seq[chr_name]
    seq = ''

    for e in exon_list:
        seq += chr_seq[e[0]:e[1] + 1]

    write_fasta(fp, gene_id, seq)

    return


def map_to_exons(old_exon, consensus_exon_list):

    for idx, cur_exon in enumerate(consensus_exon_list):
        if old_exon[0] >= cur_exon[0] and old_exon[1] <= cur_exon[1]:
            offset = old_exon[0] - cur_exon[0]

            return [cur_exon[2] + offset, cur_exon[2] + offset + (old_exon[1] - old_exon[0])], idx

    print("Warning: Can't find old_exon", old_exon, file=sys.stderr)

    return [-1, -1], -1

"""
"""
def write_transcripts_snp(fp, gene_id, trans_ids, transcripts, exon_list):

    snps = set()

    for trans_id in trans_ids:
        trans = transcripts[trans_id]
       
        if bDebug:
            pprint.pprint(trans)
        
        old_exons = trans[3]

        assert len(old_exons) >= 1
       
        new_exons = list()
        for old_exon in old_exons:
            ne, idx = map_to_exons(old_exon, exon_list)

            if bDebug:
                print(ne, idx)

            if idx == -1:
                continue

            new_exons.append(ne)


        for i in range(1, len(new_exons)):
            gap = new_exons[i][0] - new_exons[i-1][1] - 1
            if gap > 0:
                snps.add((trans[0], new_exons[i-1][1] + 1, gap, trans_id))



    snps_list = sorted(list(snps))

    snp_count = 0

    for s in snps_list:
        print('{}.{}\t{}\t{}\t{}\t{}'.format(gene_id, snp_count, 'deletion', gene_id, s[1], s[2]), file=fp)
        snp_count += 1


    return


"""
"""
def write_transcripts_map(fp, gene_id, chr_name, exon_list):
    print('>{}\t{}'.format(gene_id, chr_name), file=fp)

    for i, e in enumerate(exon_list):
        fp.write('{}:{}:{}'.format(chr_name, e[0], e[1] - e[0] + 1))
        if i % 4 == 3:
            fp.write('\n')
        else:
            fp.write('\t')

    if len(exon_list) % 4 != 0:
        fp.write('\n')

    return

"""
"""
def extract_transcript_graph(genome_file, gtf_file, base_fname):
    genome_seq = read_genome(genome_file)
    genes, transcripts = read_transcript(genome_seq, gtf_file)

    """
    if bDebug:
        pprint.pprint(genes)
        #pprint.pprint(transcripts)

        return
    """

    # get sorted id
    gene_ids = sorted(list(genes.keys()))

    with open(base_fname + ".trans.fa", "w") as trseq_file,\
        open(base_fname + ".trans.snp", "w") as trsnp_file,\
        open(base_fname + ".trans.map", "w") as trmap_file:

        # for each gene, build a consensus exons
        #   exon in consensus exons is not overlapped to other exon

        for gene_id in gene_ids:

            if bDebug:
                if gene_id != "ENSG00000244625":
                    continue

            trans_ids = genes[gene_id]

            chrom = transcripts[trans_ids[0]][0]

            tmp_exons_list = list()
            for tid in trans_ids:
                tmp_exons_list += transcripts[tid][3]

            exon_list = make_consensus_exons(tmp_exons_list)

            if bDebug:
                pprint.pprint(sorted(tmp_exons_list))
                pprint.pprint(exon_list)

            # print sequeunce
            write_transcripts_seq(trseq_file, gene_id, chrom, genome_seq, exon_list) 

            # print snps
            write_transcripts_snp(trsnp_file, gene_id, trans_ids, transcripts, exon_list)

            # print exon_map
            write_transcripts_map(trmap_file, gene_id, chrom, exon_list)



            """
            exon_count = 0
            for t_id in genes[gene_id]:
                exon_count += len(transcripts[t_id][3])
            print(gene_id, len(genes[gene_id]), exon_count)
            """


    return



if __name__ == '__main__':
    parser = ArgumentParser(
            description='Extract transcripts')

    parser.add_argument('-g', '--gtf-file',
            dest='gtf_file',
            type=FileType('r'),
            help='input GTF fie')

    parser.add_argument('-r', '--genome',
            dest='genome_file',
            type=FileType('r'),
            help='reference genome file')

    parser.add_argument('-o', '--output',
            dest='out_fname',
            type=str,
            help='output filename prefix')

    parser.add_argument('--debug',
            dest='bDebug',
            action='store_true',
            help='Run in debug mode')
    parser.add_argument('--verbose',
            dest='bVerbose',
            action='store_true',
            help='Show more messages')

    args = parser.parse_args()
    if not args.gtf_file or not args.genome_file or not args.out_fname:
        parser.print_help()
        exit(1)


    if args.bDebug != None:
        bDebug = args.bDebug

    if args.bVerbose != None:
        bVerbose = args.bVerbose


    extract_transcript_graph(
            args.genome_file,
            args.gtf_file,
            args.out_fname)

