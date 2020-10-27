HISAT-3N
============

Overview
-----------------
HISAT-3N (hierarchical indexing for spliced alignment of transcripts - 3 nucleotides)
is an ultrafast and memory-efficient sequence aligner designed for nucleotide conversion
sequencing technologies. HISAT-3N index contains two HISAT2 indexes which require memory small: 
for the human genome, it requires 9 GB for standard 3N-index and 10.5 GB for repeat 3N-index.
The repeat 3N-index could be used to align one read to thousands position 3 times faster standard 3N-index.
HISAT-3N is developed based on [HISAT2](https://github.com/DaehwanKimLab/hisat2), 
which is particularly optimized for RNA sequencing technology. 
HISAT-3N can be used for any base-converted sequencing reads include BS-seq, SLAM-seq, scBS-seq, scSLAM-seq, and TAPS.

Getting started
============
HISAT-3N requires a 64-bit computer running either Linux or Mac OS X and at least 16 GB of RAM. 

A few notes:  

1. The repeat 3N index building process requires 256 GB of RAM.
2. The standard 3N index building requires no more than 16 GB of RAM.
3. The alignment process with either standard or repeat index requires no more than 16 GB of RAM.

Install
------------
   
    git clone https://github.com/DaehwanKimLab/hisat-3n.git
    cd hisat-3n
    make

Build a 3N index with hisat-3n-build
-----------
`hisat-3n-build` builds a 3N-index, which contains two hisat2 indexes, from a set of DNA sequences. For standard 3N-index,
each index contains 16 files with suffix `.3n.*.*.ht2`.
For repeat 3N-index, there are 16 more files in addition to the standard 3N-index, and they have the suffix 
`.3n.*.rep.*.ht2`. 
These files constitute the hisat-3n index and no other file is needed to alignment reads to the reference.

* Sample argument for standard HISAT-3N index building:  
`hisat-3n-build genome.fa genome`  

* Sample argument for repeat HISAT-3N index building (require 256 GB memory):  
`hisat-3n-build --repeat-index genome.fa genome` 

It is optional to make the graph index and add SNP or spicing site information to the index, to increase the alignment accuracy.
for more detail, please check the [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/).

    # Standard HISAT-3N integrated index with SNP information
    hisat-3n-build --exons genome.exon genome.fa genome 
    
    # Standard HISAT-3N integrated index with splicing site information
    hisat-3n-build --ss genome.ss genome.fa genome 
    
    # Repeat HISAT-3N integrated index with SNP information
    hisat-3n-build --repeat-index --exons genome.exon genome.fa genome 
    
    # Repeat HISAT-3N integrated index with splicing site information
    hisat-3n-build --repeat-index --ss genome.ss genome.fa genome 

Alignment with hisat-3n
------------
After we build the HISAT-3N index, you are ready to use HISAT-3N for alignment. 
HISAT-3N uses the HISAT2 argument but has some extra arguments. Please check [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/) for more detail.

For human genome reference, HISAT-3N requires about 9GB for alignment with standard 3N-index and 10.5 GB for repeat 3N-index.

* `--base-change`  
    Provide which base is converted in the sequencing process to another base. Please enter
    2 letters separated by ',' for this argument. The first letter should be the converted base, the second letter should be
    the converted to base. For example, during slam-seq, some 'T' is converted to 'C',
    please enter `--base-change T,C`. During bisulfite-seq, some 'C' is converted to 'T', please enter `--base-change C,T`.
    If you want to align non-converted reads to the regular HISAT2 index, do not use this option.
       
* `--index/-x`  
    The index for HISAT-3N.  The basename is the name of the index files up to but not including the suffix `.3n.*.*.ht2` / etc. 
    For example, you build your index with basename 'genome' by HISAT-3N-build, please enter `--index genome`.
      
* `--repeat-limit` 
    You can set up the number of alignment will be check for each repeat alignment. You may increase the number to let hisat-3n 
    output more, if a read has multiple mapping. We suggest the repeat limit number for paired-end reads alignment is no more 
    than 1,000,000. default: 1000.

* `--unique-only` 
    Only output uniquely aligned reads.
    
Sample argument:  
* Single-end slam-seq reads (have T to C conversion) alignment with standard 3N-index:  
`--hisat-3n --index genome -f -U read.fa -S output.sam --base-change T,C`

* Paired-end bisulfite-seq reads (have T to C conversion) alignment with repeat 3N-index:   
`--hisat-3n --index genome -f -1 read_1.fa -2 read_2.fa -S output.sam --base-change C,T`

* Single-end TAB-seq reads (have T to C conversion) alignment with repeat 3N-index and only output unique aligned result:   
`--hisat-3n --index genome -q -U read.fq -S output.sam --base-change C,T --unique`

Publication
============

* HISAT-3N paper

* HIAST2 paper  
Kim, D., Paggi, J.M., Park, C. et al. [Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://www.nature.com/articles/s41587-019-0201-4). Nat Biotechnol 37, 907–915 (2019)

