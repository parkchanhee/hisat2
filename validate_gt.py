#!/usr/bin/env python3

from sys import stderr, exit
from collections import defaultdict, Counter
from argparse import ArgumentParser, FileType
import itertools
import pprint
import bisect


bDebug = False
bVerbose = False
bUniqueSNP = False
bChrTome = False


def read_snps(snp_file):
    print('Read SNPs')
    snps = defaultdict(list)

    if not snp_file:
        return snps

    for line in snp_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        try:
            snpID, type, chr, pos, data = line.split('\t')
        except ValueError:
            continue

        assert type in ["single", "deletion", "insertion"]
        if type == "deletion":
            data = int(data)

        snps[snpID] = [snpID, type, int(pos), data, chr]

    return snps


def read_haplotypes(hap_file):
    print('Read Haplotypes')

    haplotypes = defaultdict(list)

    if not hap_file:
        return haplotypes

    for line in hap_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # ht-id chrname left right snplist(,)

        try:
            htid, chrname, left, right, snplist = line.split('\t')

            haplotypes[htid] = [htid, chrname, int(left), int(right), snplist.split(',')]
        except ValueError:
            continue

    return haplotypes


def validate_haplotype(haplotype, snps):
    count = 0
    for htid, hap in haplotype.items():
        assert htid == hap[0]

        chrname = hap[1]
        left = hap[2]
        right = hap[3]
        snp_list = hap[4]

        assert snp_list

        first_snp = snps[snp_list[0]]
        last_snp = snps[snp_list[-1]]

        first_pos = first_snp[2]
        last_pos = last_snp[2]
        if last_snp[1] == 'deletion':
            last_pos = last_snp[2] + last_snp[3] - 1

        if not (left == first_pos and right == last_pos):
            print(htid, hap, first_snp, last_snp)

        count += 1

    print('Total {} hap processed'.format(count))

    return


def validate_GT(snp_file, hap_file):
    snps = read_snps(snp_file)
    haplotypes = read_haplotypes(hap_file)

    validate_haplotype(haplotypes, snps)

    return

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract transcripts')

    parser.add_argument('-s', '--snp-file',
                        dest='snp_file',
                        type=FileType('r'),
                        help='input SNP file')

    parser.add_argument('-p', '--haplotype',
                        dest='hap_file',
                        type=FileType('r'),
                        help='input Haplotype file')

    parser.add_argument('--debug',
                        dest='bDebug',
                        action='store_true',
                        help='Run in debug mode')

    parser.add_argument('--verbose',
                        dest='bVerbose',
                        action='store_true',
                        help='Show more messages')

    args = parser.parse_args()
    if not args.snp_file or not args.hap_file:
        parser.print_help()
        exit(1)

    if args.bDebug is not None:
        bDebug = args.bDebug

    if args.bVerbose is not None:
        bVerbose = args.bVerbose

    validate_GT(args.snp_file, args.hap_file)
