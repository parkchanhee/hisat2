#!/usr/bin/env python3
import sys


"""
"""
def parse_gene_line(fields, genes):
    values = fields[8]

    begin_pos = int(fields[3]) - 1
    end_pos = int(fields[4]) - 1
    strand = fields[6]
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

    gene_id = values_dict['gene_id']

    if gene_id in genes:
        print('Duplicated gene: {}'.format(gene_id))
    else:
        genes[gene_id] = [gene_id, values_dict['gene_name'], values_dict['gene_biotype'], begin_pos, end_pos, strand]

    return


def parse_transcript_line(fields, genes, transcripts):
    values = fields[8]

    values_dict = {}
    for attr in values.split(';'):
        name, _, value = attr.strip().partition(' ')
        values_dict[name] = value.strip('"')

    gene_id = values_dict['gene_id']
    transcript_id = values_dict['transcript_id']

    if gene_id not in genes:
        print('Wrong transcript: {} / {}'.format(gene_id, transcript_id))
        return

    if transcript_id in transcripts:
        print('Duplicated transcript: {}  {}'.format(gene_id, transcript_id))
        return

    gene_biotype = values_dict['gene_biotype']
    transcript_biotype = values_dict['transcript_biotype']
    transcripts[transcript_id] = [list(), gene_biotype, transcript_biotype]

    return

gtf_file = sys.argv[1]

genes = {}
transcripts = {}

with open(gtf_file) as fp:
    for line in fp:
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        try:
            fields = line.split('\t')

        except ValueError:
            print('Wrong line:', line)
            continue

        feature = fields[2]
        if feature == "gene":
            parse_gene_line(fields, genes)
        elif feature == 'transcript':
            parse_transcript_line(fields, genes, transcripts)
        else:
            continue

#print(genes)

def genes_groupby_biotype(genes):
    group_biotype = {}

    for key, value in genes.items():
        biotype = value[2]

        if biotype in group_biotype:
            group_biotype[biotype].append(key)
        else:
            group_biotype[biotype] = [key]

    return group_biotype


print('Total genes: {}'.format(len(genes)))

print('Number of genes by gene_biotype')
groupby_biotype = genes_groupby_biotype(genes)
for key, value in groupby_biotype.items():
    print('\t{}: {}'.format(key, len(value)))

v = sum([len(value) for value in groupby_biotype.values()])
assert v == len(genes)

print(len(transcripts))


for gene in genes.items():
    print(gene)


print('-----------------')

sorted_gene = sorted(genes, key=lambda gene: (gene[3], gene[4], gene[5])) # sort by start_pos, end, strand

for gene in genes.items():
    print(gene)


"""
for gene in genes:
    print(gene[2])
"""