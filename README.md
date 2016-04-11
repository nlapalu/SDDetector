# SDDetector
Segmental Duplication Detection tool

## Author

## Install

### Install SSDetector

prerequesite: 

* Python 2.7.X
* BioPython, only if you want to parse Blast results from xml format

### Additional tools

* bedtools
* ncbiblast+


## Example

genome.fasta is your genome in multi-fasta file
genome_TE.gff is the Transposable, Repeat elements annotation file in gff3  

### Build a masked blast database

__Convert the genome fasta file in a soft-masked fasta file (upper cases to lower cases)__

maskFastaFromBed -fi genome.fasta -fo genome_masked.fasta -bed genome_TE.gff -soft

__Index your genome with masking information__

convert2blastmask -in genome_masked.fasta -parse_seqids -masking_algorithm REPET -masking_options "REPET, URGI" -outfmt maskinfo_asn1_bin -out genome_masked.asnb
makeblastdb -dbtype nucl -in genome_masked.fasta -out genome_masked -parse_seqids -mask_data genome_masked.asnb

__Check your masking info:__

blastdbcmd -db genome_masked -info

Database: /tmp/genome_masked.fasta
    28 sequences; 50,819,261 total bases

Date: Feb 23, 2016  4:05 PM    Longest sequence: 6,042,495 bases

Available filtering algorithms applied to database sequences:

Algorithm ID  Algorithm name      Algorithm options                       
100           other               REPET, URGI                             

Volumes:
    /tmp/genome_masked


### Perfom blast analysis

__Blast in XML format__

blastn -num_threads 2 -task megablast -db seq -query seq.fasta -out blast.xml -outfmt 5 -db_soft_mask 100

__Blast in tab-delimited format with required fields__

blastn -num_threads 2 -task megablast -db seq -query seq.fasta -out blast.tab -outfmt "6 qseqid sseqid qstart qend sstart send length nident" -db_soft_mask 100

### Detect Segmental Duplications

## How to cite

If you use SDDetector, please cite:

## References


