#!/usr/bin/env python

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
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

from __future__ import print_function

from sys import stderr, exit
from collections import defaultdict as dd, Counter
from argparse import ArgumentParser, FileType
import bisect

# global
flagLoadGenome=True
verbose=False
debug=False

class Range:
    def __init__(self, left, right):
        self.left = left
        self.right = right
        return

    def __repr__(self):
        return '({}, {})'.format(self.left, self.right)


def intersect(a, b):
    left = a
    right = b

    if b[0] < a[0]:
        left = b
        right = a

    if left[1] < right[0]:
        return None

    return (right[0], min(left[1], right[1]))

def union(a, b):
    assert intersect(a, b) != None
    return (min(a[0], b[0]), max(a[1], b[1]))


class NonOverlappingRangeSet:
    def __init__(self):
        self.ranges = []
        return

    def put(self, point):
        # point is all-inclusive  [a, b]
        assert point[0] <= point[1]

        if len(self.ranges) == 0:
            self.ranges.append(point[0])
            self.ranges.append(point[1])
        else:
            # check overlapping range
            overlap_list = []

            left_idx = bisect.bisect_left(self.ranges, point[0])
            right_idx = bisect.bisect_left(self.ranges, point[1])

            #print left_idx, right_idx

            if left_idx == 0 and right_idx == 0:
                if point[1] < self.ranges[right_idx]:
                    self.ranges.insert(0, point[1])
                    self.ranges.insert(0, point[0])
                    return
                else:
                    # point[1] == ranges[right_idx]
                    # merge point with (ranges[0], ranges[1])
                    
                    assert point[0] <= self.ranges[0]
                    assert point[0] <= self.ranges[1]
                    self.ranges[0] = point[0]

                    return

            if left_idx == len(self.ranges):
                # append point to ranges
                self.ranges.append(point[0])
                self.ranges.append(point[1])
                return

            if right_idx == len(self.ranges):
                # range[l - 1] < point[0] <= range[l] ... ranges[-1] < point[1] 
                # if l is an even
                #   new_left is l, ranges[new_left] = point[0]
                # if l is an odd
                #   new_left is l-1, values aren't changes
                #
                # right of new ranges is len(ranges)-1, ranges[-1] = point[1]
                #
                if left_idx % 2 == 0:
                    new_left = left_idx
                    self.ranges[new_left] = point[0]
                else:
                    assert left_idx > 0
                    new_left = left_idx - 1

                new_right = right_idx - 1
                self.ranges[new_right] = point[1]

                assert new_left % 2 == 0
                assert new_right % 2 == 1

                tmp_ranges = self.ranges[0:new_left+1] + self.ranges[new_right:]
                self.ranges = tmp_ranges

                return



            if left_idx == right_idx:
                # m = left_idx = right_idx > 0 
                # -> range[m-1] < point[0] <= point[1] <= range[m]
                # if m is an even, 
                #    if point[1] < range[m]
                #       insert new point at m
                #    if point[1] == range[m]
                #       point and (range[m], range[m+1]) can be merged
                #           -> point[m] = point[0] 
                # if m is an odd,
                #    (range[m-1], range[m]) merge point
                m = left_idx
                if m % 2 == 0:
                    # even
                    if point[1] < self.ranges[m]:
                        self.ranges.insert(m, point[1])
                        self.ranges.insert(m, point[0])
                    else:
                        self.ranges[m] = point[0]
                else:
                    # odd
                    pass
                return

            assert left_idx < right_idx
            
            # l = left_idx
            # r = right_idx

            # ranges[l - 1] < point[0] <= ranges[l] <= ranges[r - 1] < point[1] <= ranges[r]
            #
            # if l is an even
            #   left of new range is l, ranges[l] = point[0]
            # if l is an odd
            #   left of new range is l-1, values aren't change
            #
            # if r is an odd
            #   right of new range is r, ranges[r] are not changed
            # if r is an even
            #   if point[1] == ranges[r]
            #     right of new range is ranges[r+1], values aren't changed
            #   if point[1] < ranges[r]
            #     right of new range is ranges[r-1] = point[1]
            #           
            # delete new_left + 1 ~ new_right - 1
            #

            new_left = left_idx
            new_right = right_idx

            if left_idx % 2 == 0:
                new_left = left_idx
                self.ranges[new_left] = point[0]
            else:
                new_left = left_idx - 1

            if right_idx % 2 == 0:
                if self.ranges[right_idx] == point[1]:
                    new_right = right_idx + 1
                    assert new_right < len(self.ranges)
                else:
                    new_right = right_idx - 1
                    self.ranges[new_right] = point[1]
            else:
                new_right = right_idx

            # delete
        
            assert new_left < new_right
            #assert (new_left + 1 ) <= (new_right - 1)

            tmp_ranges = self.ranges[0:new_left+1] + self.ranges[new_right:]

            self.ranges = tmp_ranges

            assert len(self.ranges) % 2 == 0
        return

    def __repr__(self):
        return str(self.ranges)

    def __len__(self):
        return len(self.ranges)

    def show_ranges(self):
        for pt in get_ranges(self):
            print('rs', pt[0], pt[1])
        return

    def get_ranges(self):
        return [list(a) for a in zip(self.ranges[::2], self.ranges[1::2])]


def write_fasta(fasta_file, ref_id, seq, width=60):
    print('>{}'.format(ref_id), file=fasta_file)

    length=len(seq)
    p=0

    while p < length:
        print(seq[p:p+width], file=fasta_file)
        p += width

"""
"""
def read_genome(genome_file):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence

    return chr_dic



"""
"""
def write_transcripts_seq(trseq_file, chr_dic, chrom, exon_maps):
    ref_id = chrom + '_transcriptome'
    seq = ''

    if chrom not in chr_dic:
        return

    for (offset, (left, right)) in exon_maps:
        seq += chr_dic[chrom][left:right+1]

    write_fasta(trseq_file, ref_id, seq)


def write_transcripts_snp(trsnp_file, genes, trans, exon_maps):

    snps = set()
    for tran_key, transcripts in trans.items():
        chrom = tran_key[0]

        new_exons = map_exons(transcripts[2], exon_maps[chrom])

        for i in range(len(new_exons) - 1):
            ex = new_exons[i]
            nex = new_exons[i + 1]

            if nex[0] - ex[1] == 1:
                continue

            deletion_start = ex[1] + 1
            deletion_size = nex[0] - ex[1] - 1 

            assert deletion_size > 0

            snps.add((chrom, deletion_start, deletion_size))



    """
    assert len(tran_keys) > 0

    snps = set()
    for tk in tran_keys:
        exons = trans[tk][2]
        new_exons = map_exons(exons, exon_maps)

        for i in range(len(new_exons) - 1):
            ex = new_exons[i]
            nex = new_exons[i + 1]

            if nex[0] - ex[1] == 1:
                continue

            deletion_start = ex[1] + 1
            deletion_size = nex[0] - ex[1] - 1 

            assert deletion_size > 0

            snps.add((deletion_start, deletion_size))
    """

    # sort by pos
    snps = list(snps)
    snps.sort()
    # print
    idx=0
    for snp in snps:
        snp_id = 'tran.' + str(idx)
        idx += 1
        print("{}\t{}\t{}\t{}\t{}".format(snp_id, "deletion", snp[0] + '_transcriptome', snp[1], snp[2]), file=trsnp_file)


def write_transcripts_map(trmap_file, genes, trans, chrom, exon_map):

    item_per_line = 4
    num_items = len(exon_map)

    last_exon = exon_map[-1]
    total_size = last_exon[0] + (last_exon[1][1] - last_exon[1][0] + 1)

    print(">{}\t{}\t{}".format(chrom + '_transcriptome', num_items, total_size), file=trmap_file)

    #chrom = tran_keys[0][0]
    i = 0

    #print(exon_maps)

    while i < num_items:
        buf = ""
        j = i
        while j < min(num_items, i + item_per_line):
            off, exon  = exon_map[j]
            size = exon[1] - exon[0] + 1
            buf += "{}:{}:{}\t".format(chrom, exon[0], size)
            j += 1

        # trim last \t
        buf = buf[:-1]

        print(buf, file=trmap_file)
        i += item_per_line


def merge_exons(gene_id, tran_keys, trans):

    gene_exon_set = set()
    t_chrom = None
    t_strand = None

    exon_count = 0
    for tran_key in tran_keys:
        chrom, strand, exons = trans[tran_key]
        transcript_id = tran_key[1]

        if t_chrom is None:
            t_chrom = chrom
            t_strand = strand
        else:
            assert t_chrom == chrom
            assert t_strand == strand

        for exon in exons:
            ekey = (exon)
            exon_count += 1
            gene_exon_set.add(ekey)

    if verbose:                
        print("Gene ID: {}, Total Exon: {}, Uniq Exon: {}, Total Trascripts: {}".format(gene_id, exon_count, len(gene_exon_set), len(tran_keys)), file=stderr)

    gene_exon_list = list(gene_exon_set)
    gene_exon_list.sort()
    #print(gene_exon_list)

    tran_exon_index = {}

    # rebuild transcript
    for tran_key in tran_keys:
        chrom, strans, exons = trans[tran_key]
        exon_id_list = [gene_exon_list.index(exon) for exon in exons]
        if debug:
            print(tran_key[1], exon_id_list) 

        tran_exon_index[tran_key[1]] = exon_id_list
                
        for i in range(len(exon_id_list)-1):
            assert exon_id_list[i] < exon_id_list[i+1]

    return gene_exon_list, tran_exon_index 

"""
"""
def write_transcriptome_fasta(genes, trans, chr_dic, trseq_file):
    for gene_id, tran_keys in genes.items():
        # tran_key = (chrom, transcript_id)
        for tran_key in tran_keys:
            chrom, strand, exons = trans[tran_key]
            seq = ''
            ref_id = tran_key[1] + '^' + chrom
            if chrom not in chr_dic:
                continue

            for exon in exons:
                left=exon[0]-1
                right=exon[1]-1
                seq += chr_dic[chrom][left:right+1]

                #print('{}\t{}\t{}\t{}\t{}'.format(
                #            gene_id, tran_key[1], exon[0], exon[1], seq), file=stderr)
                
            write_fasta(trseq_file, ref_id, seq) 

def write_transcriptome_map(genes, trans, chr_dic, trmap_file):
    return


def build_consensus_exons(gene_id, tran_keys, trans):

    rs = NonOverlappingRangeSet()
    for tk in tran_keys:
        # tk = [chrom, transcript_id]
        
        # [chrom, strand, exons]
        transcript = trans[tk] 
        for ex in transcript[2]:
            if debug:
                print(tk[1][-6:], ex[0], ex[1])
            rs.put(ex)

    tmp_exons = rs.get_ranges()
    exon_maps = []
    offset = 0
    for te in tmp_exons:
        exon_maps.append([offset, te])
        offset += te[1] - te[0] + 1

    if debug:
        rs.show_ranges()

    return exon_maps 

def show_exons(gene_id, tran_keys, trans):

    exons_counter = Counter()
    #tran_keys.sort()
    for tk in tran_keys:
        # tk = [chrom, transcript_id]
        
        # [chrom, strand, exons]
        transcript = trans[tk] 
        for ex in transcript[2]:
            #print(ex)
            exons_counter[ex] += 1

    #print(exons_counter)
    #print(exons_counter.values())
    print(gene_id, 'Total transcripts:', len(tran_keys))
    print(gene_id, 'Total exons:', sum(exons_counter.values()))
    print(gene_id, 'Unique exons:', len(exons_counter))

    # merging exons
    #rs = NonOverlappingRangeSet()
    #rs.put((1,10))
    #rs.put((20,30))
    #rs.put((14,19))
    #rs.put((18,21))
    #print(rs)

    rs = NonOverlappingRangeSet()

    for tk in tran_keys:
        # tk = [chrom, transcript_id]
        
        # [chrom, strand, exons]
        transcript = trans[tk] 
        for ex in transcript[2]:
            if debug:
                print(tk[1][-6:], ex[0], ex[1])
            rs.put(ex)

    #print(rs)
    if debug:
        rs.show_ranges()
    print(gene_id, 'Merged exons:', len(rs)) 


    tmp_exons = rs.get_ranges()
    exon_maps = []
    offset = 0
    for te in tmp_exons:
        exon_maps.append([offset, te])
        offset += te[1] - te[0] + 1

    print(exon_maps)

    # map old exons to new sequence
    for tk in tran_keys:
        # tk = [chrom, transcript_id]
        
        exons = trans[tk][2]
        # [chrom, strand, exons]
        new_exons = map_exons(exons, exon_maps) 
        print('Old:', exons)
        print('New:', new_exons)



    """
    print(intersect((10, 20), (30, 40))) # None
    print(intersect((30, 40), (10, 20))) # None
    print(intersect((10, 20), (10, 20))) # (10, 20)
    print(intersect((10, 20), (15, 30))) # (15, 20)
    print(intersect((10, 20), (1, 15))) # (10, 15)
    print(intersect((10, 20), (11, 15))) # (11, 15)
    print(intersect((10, 20), (1, 30))) # (10, 20)
    """
    return


def make_common_exons(genes, trans, chr_dic, fp_seq, fp_map, fp_snp):

    # 
    rs_dic = {}

    # for each gene
    # 1. build consensus exons
    # 2. get new exon position
    # 3. write files
    for gene_id, tran_keys in genes.items():
        if debug:
            if gene_id != "ENSG00000244625":
                continue
        
        #show_exons(gene_id, tran_keys, trans) 
        #exon_maps = build_consensus_exons(gene_id, tran_keys, trans)
        #print(gene_id, exon_maps)
        
        for tk in tran_keys:
            chrom, tran_id = tk

            if chrom not in rs_dic:
                rs_dic[chrom] = NonOverlappingRangeSet()

            transcripts = trans[tk]
            for exon in transcripts[2]:
                if debug:
                    print(chrom + ':' + tk[1][-6:], exon[0], exon[1])
                rs_dic[chrom].put(exon)

    exon_maps = {}

    for chrom, rs in rs_dic.items():
        tmp_exons = rs.get_ranges()

        exon_maps[chrom] = []
        offset = 0
        for te in tmp_exons:
            exon_maps[chrom].append([offset, te])
            offset += te[1] - te[0] + 1

        if debug:
            rs.show_ranges()

    for chrom, exon_map in exon_maps.items():
        # write consensus exons
        write_transcripts_seq(fp_seq, chr_dic, chrom, exon_map)

    # for each transcript, generate snps
    write_transcripts_snp(fp_snp, genes, trans, exon_maps)

    # write exon map
    for chrom, exon_map in exon_maps.items():
        write_transcripts_map(fp_map, genes, trans, chrom, exon_map) 

    return


def map_exons(exons, exon_maps):

    new_exons = []

    for ex in exons:
        left, right = ex
        length = right - left + 1

        found = False
        for m in exon_maps:
            exon_range = m[1]
            if left >= exon_range[0] and right <= exon_range[1]:
                # found
                found = True
                break
        assert found


        exon_offset, exon_range = m
        offset = left - exon_range[0]
        new_left = exon_offset + offset
        new_right = new_left + length - 1

        #print(ex, (new_left, new_right))

        new_exons.append((new_left, new_right))

    return new_exons 


def extract_transcriptome(gtf_file, ref_file, out_name):
    trseq_file = open(out_name + ".transcriptome.fa", "w")
    trsnp_file = open(out_name + ".transcriptome.snp", "w")
    trmap_file = open(out_name + ".transcriptome.map", "w")

    genes = dd(list)
    trans = {}

    if flagLoadGenome:
        chr_dic = read_genome(ref_file)
    else:
        chr_dic = {}

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            continue
        left, right = int(left) - 1, int(right) - 1

        #print(chrom, feature, left, right)
        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';')[:-1]:
            attr, _, val = attr.strip().partition(' ')
            values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            continue

        transcript_id = values_dict['transcript_id']
        gene_id = values_dict['gene_id']
        key = (chrom, transcript_id)
        if key not in trans:
            trans[key] = [chrom, strand, [(left, right)]]
            genes[gene_id].append(key)
        else:
            trans[key][2].append((left, right))
            assert chrom == trans[key][0]
            assert strand == trans[key][1]



    merged_count = 0
    # Sort exons and merge where separating introns are <=5 bps
    for (chrom, tran), [chrom, strand, exons] in trans.items():
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                assert exons[i][0] > exons[i-1][1]

                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tm = tmp_exons[-1]
                    del tmp_exons[-1]
                    tmp_exons.append((tm[0], exons[i][1]))
                    #tmp_exons[-1][1] = exons[i][1]
                    merged_count += 1
                else:
                    tmp_exons.append(exons[i])
            trans[(chrom, tran)] = [chrom, strand, tmp_exons]

    print('Merged Count:', merged_count, file=stderr)


    # Write genes, trans
    #write_genes(genes, trans, chr_dic, trseq_file, trmap_file, trsnp_file)

    # Write transcriptome fasta file
    #write_transcriptome_fasta(genes, trans, chr_dic, trseq_file)

    # Write transcriptome map file
    #write_transcriptome_map(genes, trans, chr_dic, trmap_file)

    make_common_exons(genes, trans, chr_dic, trseq_file, trmap_file, trsnp_file)

    if verbose:
        exon_lengths, intron_lengths, trans_lengths = \
            Counter(), Counter(), Counter()
        for chrom, strand, exons in trans.values():
            tran_len = 0
            for i, exon in enumerate(exons):
                exon_len = exon[1]-exon[0]+1
                exon_lengths[exon_len] += 1
                tran_len += exon_len
                if i == 0:
                    continue
                intron_lengths[exon[0] - exons[i-1][1]] += 1
            trans_lengths[tran_len] += 1

        print('genes: {}, genes with multiple isoforms: {}'.format(
                len(genes), sum(len(v) > 1 for v in genes.values())),
              file=stderr)
        print('transcripts: {}, transcript avg. length: {:d}'.format(
                len(trans), sum(trans_lengths.elements())/len(trans)),
              file=stderr)
        print('exons: {}, exon avg. length: {:d}'.format(
                sum(exon_lengths.values()),
                sum(exon_lengths.elements())/sum(exon_lengths.values())),
              file=stderr)
        print('introns: {}, intron avg. length: {:d}'.format(
                sum(intron_lengths.values()),
                sum(intron_lengths.elements())/sum(intron_lengths.values())),
              file=stderr)
        print('average number of exons per transcript: {:d}'.format(
                sum(exon_lengths.values())/len(trans)),
              file=stderr)

    trseq_file.close()
    trmap_file.close()
    trsnp_file.close()

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract transcriptome')
    parser.add_argument('-g', '--gtf-file',
        dest='gtf_file',
        type=FileType('r'),
        help='input GTF file')
    parser.add_argument('-r', '--reference',
        dest='ref_file',
        type=FileType('r'),
        help='reference genome file')
    parser.add_argument('-o', '--output',
        dest='out_name',
        type=str,
        help='output filename prefix')
    """
    parser.add_argument('-t', '--transcriptome',
        dest='trseq_file',
        type=FileType('w'),
        help='[OUT] transcriptome sequences file(fasta)')
    parser.add_argument('-m', '--tr-map',
        dest='trmap_file',
        type=FileType('w'),
        help='[OUT] transcriptome mapping file')
    """
    parser.add_argument('--without-genome',
        dest='flagLoadGenome',
        action='store_false',
        help="Don't loading genome reference")
    parser.add_argument('--debug',
        dest='flagDebug',
        action='store_true',
        help="Run in debug mode")
    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='also print some statistics to stderr')

    args = parser.parse_args()
    if not args.gtf_file or not args.ref_file or not args.out_name:
        parser.print_help()
        exit(1)

    if args.flagLoadGenome != None:
        flagLoadGenome = args.flagLoadGenome
    if args.flagDebug != None:
        debug = args.flagDebug

    if args.verbose != None:
        verbose = args.verbose

    extract_transcriptome(
            args.gtf_file,
            args.ref_file,
            args.out_name)
