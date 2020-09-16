HISAT-3N is under development!
============

What is HISAT-3N?
-----------------

HISAT-3N is a three nucleotide aligner integrated in HISAT2. 
It could be used for any base-converted sequencing reads like
slam-seq, bisulfite-seq.


Build index for HISAT-3N
============
hisat2-build
-----------
You need to run hisat-3n-build to build HISAT-3N index.

Sample argument for HISAT-3N index building without repeat index:  
`hisat-3n-build -p 10 genome.fa genome`  

Sample argument for HISAT-3N index building with repeat index:  
`hisat-3n-build -p 10 --auto-repeat-index 100-300 genome.fa genome` 

HISAT-3N index with splicing site information is under development!  
PLEASE DO NOT USE THIS!

    HISAT2_SS_SCRIPT=./hisat2_extract_splicesites.py
    HISAT2_EXON_SCRIPT=./hisat2_extract_exons.py
    download ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
    
    tar zxvf Homo_sapiens.GRCh38.84.gtf.gz
    mv  Homo_sapiens.GRCh38.84.gtf genome.gtf
    GTF_FILE=genome.gtf
    ${HISAT2_SS_SCRIPT} ${GTF_FILE} > genome.ss
    ${HISAT2_EXON_SCRIPT} ${GTF_FILE} > genome.exon

Sample argument for HISAT-3N index building with splicing site information:  
`hisat-3n-build genome.fa genome -p 10 --ss genome.ss --exons genome.exon`   
 


Use HISAT-3N for alignment
============
After we build HISAT-3N index, you are ready to use HISAT-3N for alignemnt.
HISAT-3N use HISAT2 argument but has some extra argument.

For human genome reference, HISAT-3N took about 9GB for alignment.

* `--base-change`  
    Provide which base is converted in sequencing process to other base. Please enter
    2 letter seperated by ',' for this argument. First letter should be the converted base, second letter should be
    the converted to base. For example, during slam-seq, some 'T' is converted to 'C',
    please enter `--base-change T,C`. During bisulfite-seq, some 'C' is converted to 'T', please enter `--base-change C,T`.
 
* `--no-base-change`    
    Align reads as regular hisat2. 
       
* `--index/-x`  
    The index for HISAT-3N.  The basename is the name of the index files up to but not including the final `TLA.1.1.ht2` / etc. 
    For example, you build your index with basename 'genome' by HISAT-3N-build, please enter `--index genome`.
      
* `--repeat-limit` 
    You could set up number of alignment will be check for each repeat alignment. default: 1000.

Sample argument:  
* If you want to align your slam-seq reads:  
`--hisat-3n  --index genome -f -1 read_1.fa -2 read_2.fa -S output.sam --base-change T,C -p 10`

* If you want to align your bisulfite-seq reads:   
`--hisat-3n  --index genome -f -1 read_1.fa -2 read_2.fa -S output.sam --base-change C,T -p 10`

