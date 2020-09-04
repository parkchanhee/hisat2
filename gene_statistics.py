#!/usr/bin/env python3
import sys


bVerbose = False

coding_geneset = [
    "protein_coding",
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_V_gene",
    "polymorphic_pseudogene", # ??
]

non_gene_geneset = [
    "TEC"
]

misc_non_coding_geneset = [
    "misc_RNA",
    "ribozyme",
    "Mt_tRNA",
    "Mt_rRNA",

]

small_non_coding_geneset = [
    "snRNA",
    "rRNA",
    "miRNA",
    "scaRNA",
    "snoRNA",
    "vault_RNA",
    "scRNA",
    "sRNA",

]

long_non_coding_geneset = [
    "lincRNA",
    "antisense",
    "sense_intronic",
    "sense_overlapping",
    "processed_transcript",
    "lncRNA",
]

pseudogene_geneset = [
    "unprocessed_pseudogene",
    "processed_pseudogene",
    "unitary_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "transcribed_processed_pseudogene",
    "transcribed_unitary_pseudogene",
    "translated_processed_pseudogene",
    "translated_unprocessed_pseudogene",
    "rRNA_pseudogene",
    "IG_C_pseudogene",
    "IG_D_pseudogene",
    "IG_J_pseudogene",
    "IG_V_pseudogene",
    "IG_pseudogene",
    "TR_J_pseudogene",
    "TR_V_pseudogene",
    "pseudogene",
]

"""
"""
def read_genome(genome_file, chr_filter = None):
    chr_dic = {}
    #chr_filter = [str(x) for x in list(range(1, 23)) + ['X', 'Y']]
    #chr_filter = None

    chr_name, sequence = "", ""
    for line in genome_file:
        if line[0] == ">":
            if chr_name and sequence and (not chr_filter or chr_name in chr_filter):
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line[:-1]

    if chr_name and sequence and (not chr_filter or chr_name in chr_filter):
        chr_dic[chr_name] = sequence

    if bVerbose:
        lengths = [len(value) for value in chr_dic.values()]
        print('Number of Chromosomes: {}'.format(len(chr_dic)))
        print('Total length: {}'.format(sum(lengths)))

    return chr_dic

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
    transcripts[transcript_id] = [list(), gene_biotype, transcript_biotype, gene_id]

    return

gtf_file = sys.argv[1]
ref_file = sys.argv[2]

print('Read genomes')
with open(ref_file) as fp:
    ref_genome = read_genome(fp)


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

def count_genes_by_group(genes_groupby_biotype, category):
    count = 0
    for cat in category:
        if cat in genes_groupby_biotype:
            #print('{}:{}'.format(cat, len(genes_groupby_biotype[cat])))
            count += len(genes_groupby_biotype[cat])
        #else:
        #    print('{} is not exist'.format(cat))

    return count


"""
"""
def merge_range(A, B):
    assert A[0] <= B[0]

    if (A[1]+1) < B[0]:
        return None

    return [A[0], max(A[1], B[1])]


def make_unoverlapped_range(ranges):
    ranges = sorted(ranges)

    cur_range = ranges[0]
    assert cur_range[0] <= cur_range[1]
    newrange = [cur_range]

    for i in range(1, len(ranges)):
        last_range = newrange[-1]
        cur_range = ranges[i]

        m = merge_range(last_range, cur_range)

        if m is None:
            newrange.append([*cur_range])
        else:
            newrange[-1] = m

    return newrange

print('Total genes in file: {}'.format(len(genes)))

print('Number of genes by gene_biotype')
groupby_biotype = genes_groupby_biotype(genes)

for key, value in groupby_biotype.items():
    print('\t{}: {}'.format(key, len(value)))

v = sum([len(value) for value in groupby_biotype.values()])
assert v == len(genes)


print('Number of genes')
print('\tCoding genes: {}'.format(count_genes_by_group(groupby_biotype, coding_geneset)))
small_non = count_genes_by_group(groupby_biotype, small_non_coding_geneset)
long_non = count_genes_by_group(groupby_biotype, long_non_coding_geneset)
misc_non = count_genes_by_group(groupby_biotype, misc_non_coding_geneset)
print('\tNon coding genes: {}'.format(small_non + long_non + misc_non))
print('\t\tSmall non coding genes: {}'.format(small_non))
print('\t\tLong non coding genes: {}'.format(long_non))
print('\t\tMisc non coding Genes: {}'.format(misc_non))
print('\tPseudogenes: {}'.format(count_genes_by_group(groupby_biotype, pseudogene_geneset)))
print('\tNon genes: {}'.format(count_genes_by_group(groupby_biotype, non_gene_geneset)))

print('Number of transcripts: {}'.format(len(transcripts)))


"""
for gene in genes.items():
    print(gene)
"""

"""
# Show first item
print(next(iter(genes.items())))
print(len(genes))
print('-----------------')
"""
sorted_gene = sorted(genes.items(), key=lambda gene: (gene[1][3], gene[1][4], gene[1][5])) # sort by start_pos, end, strand

# print(sorted_gene[0])
# print(sorted_gene[1])
# print(sorted_gene[2])
# print(len(sorted_gene))


# find regions of genes

regions = list()
for gene_id, values in genes.items():
    left_pos, right_pos = values[3:5]
    regions.append([left_pos, right_pos])

#regions.sort()

sum_all_gene = 0

for region in regions:
    sum_all_gene += (region[1] - region[0] + 1)


unoverlapped_regions = make_unoverlapped_range(regions)
sum_all_unoverlapped = 0
for region in unoverlapped_regions:
    sum_all_unoverlapped += (region[1] - region[0] + 1)


# get coding genes
coding_genes_regions = list()
for gene_id, values in genes.items():
    gene_biotype = values[2]
    left_pos, right_pos = values[3:5]
    if gene_biotype in coding_geneset:
        coding_genes_regions.append([left_pos, right_pos])

coding_genes_regions = make_unoverlapped_range(coding_genes_regions)
sum_all_coding_gene_unoverlapped = 0
for region in coding_genes_regions:
    sum_all_coding_gene_unoverlapped += (region[1] - region[0] + 1)



chr_len = len(ref_genome['22'])

print('Genome Length: {}'.format(chr_len))
print('Sum of gene regions: {}'.format(sum_all_gene))
print('Sum of gene regions(unoverlapped): {}'.format(sum_all_unoverlapped))
print('Sum of coding gene regions(unoverlapped): {}'.format(sum_all_coding_gene_unoverlapped))

# Get unique gene names
"""
unique_gene = set()
for item in sorted_gene:
    unique_gene.add(item[1][1])
print('Total items in set:', len(unique_gene))

dup_gene_list = list()
dup_genes = dict()

for item in sorted_gene:
    gene_name = item[1][1]
    gene_id = item[0]
    if gene_name in dup_genes:
        dup_genes[gene_name][0] += 1
        dup_genes[gene_name][1].append(gene_id)
    else:
        dup_genes[gene_name] = [1, [gene_id]]

c = 0
for name, values in dup_genes.items():
    if values[0] > 1:
        dup_gene_list.append(tuple([name, values[1]]))
        c += values[0]

print('Dup gene count:', len(dup_gene_list), c)
print(dup_gene_list)


# for each dup genes
for item in dup_gene_list:
    gene_name = item[0]
    gene_ids = item[1]
    print('Gene name:', gene_name)
    for gene_id in gene_ids:
        print('\t', genes[gene_id])
"""

"""
for gene in genes:
    print(gene[2])
"""

"""
for transcript in transcripts.items():
    gene_biotype = transcript[1][1]
    transcript_biotype = transcript[1][2]

    if gene_biotype == 'protein_coding' or \
        gene_biotype == 'lncRNA':
        continue

    if gene_biotype != transcript_biotype:
        print(gene_biotype, transcript_biotype)
"""
