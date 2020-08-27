#!/usr/bin/env python3
import sys


gtf_file = sys.argv[1]

genes = {}

with open(gtf_file) as fp:
    for line in fp:
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')

        except ValueError:
            print('Wrong line:', line)
            continue

        values_dict = {}
        for attr in values.split(';'):
            name, _, value = attr.strip().partition(' ')
            values_dict[name] = value.strip('"')

        """
        print('{}, {}, {}'.format(
                    values_dict['gene_id'],
                    values_dict['gene_name'],
                    values_dict['gene_biotype']))
        """

        genes.append([
                    values_dict['gene_id'],
                    values_dict['gene_name'],
                    values_dict['gene_biotype']
                    ])

#print(genes)

# for set
tmp_gene = set()
for gene in genes:
    if gene[2] == 'protein_coding':
        tmp_gene.add(gene[1])


print(len(tmp_gene))

"""
for gene in genes:
    print(gene[2])
"""
