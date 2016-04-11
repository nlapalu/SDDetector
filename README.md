# SDDetector: a segmental duplication detection tool

SDDetector has been developed to detect segmental duplications in complete genomes. The principle is based on the bioinformatic protocol proposed by *Kahja et al.* Segmental duplications are defined as regions having a sequence similarity greater than 90% and a length greater than 5000 nt. Regions could be fragmented due to inversion/insertion/deletion events, so a maximium gap of 3000 nt is allowed between fragments. Then, fragments are chained together and the chains with required criteria are reported as portential duplication regions. To fit with your genome specificity, SDDetector allows parameters modification to increase/reduce the sequence similarity threshold, the minimal length of the regions and the maximal gap size.
For an efficient detection, transposable and repetitive elements must be masked before sequence similarity search. If you do not provide a repetitive element annotation, results will contain a lot of false-positive regions. We recommend to perform a TE detection with the REPET package(*Flutre et al.*) or at least a minimal masking step with RepetMasker(*Smit et al.*).

SSDetector is developed by Nicolas Lapalu at [INRA-BIOGER](http://www.versailles-grignon.inra.fr/bioger). Please do not hesitate to contact me (nlapalu at versailles dot inra dot fr) if you have any comments or questions.

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

`maskFastaFromBed -fi genome.fasta -fo genome_masked.fasta -bed genome_TE.gff -soft`

__Index your genome with masking information__

`convert2blastmask -in genome_masked.fasta -parse_seqids -masking_algorithm REPET -masking_options "REPET, URGI" -outfmt maskinfo_asn1_bin -out genome_masked.asnb`

`makeblastdb -dbtype nucl -in genome_masked.fasta -out genome_masked -parse_seqids -mask_data genome_masked.asnb`

__Check your masking info:__

`blastdbcmd -db genome_masked -info`


```
Database: /tmp/genome_masked.fasta
    28 sequences; 50,819,261 total bases

Date: Feb 23, 2016  4:05 PM    Longest sequence: 6,042,495 bases

Available filtering algorithms applied to database sequences:

Algorithm ID  Algorithm name      Algorithm options                       
100           other               REPET, URGI                             

Volumes:
    /tmp/genome_masked
```

### Perfom blast analysis

__Blast in XML format:__

`blastn -num_threads 2 -task megablast -db genome_masked -query genome.fasta -out blast.xml -outfmt 5 -db_soft_mask 100`

__Blast in tab-delimited format with required fields:__

`blastn -num_threads 2 -task megablast -db genome_masked -query genome.fasta -out blast.tab -outfmt "6 qseqid sseqid qstart qend sstart send length nident" -db_soft_mask 100`

### Detect Segmental Duplications

__Detection from xml format:__

`segmental_duplication_detector.py blast.xml xml sdd_0.9_3000_5000.gff3 :memory: -g 3000 -l 5000 -a`

__Detection from tab-delimited format:__:

`segmental_duplication_detector.py blast.tab tab sdd_0.9_3000_5000.gff3 :memory: -g 3000 -l 5000 -a`

## How to cite

If you use SDDetector, please cite:

*__SDDetector, Lapalu.N, https://github.com/nlapalu/SDDetector__*

## References

* Flutre T, Duprat E, Feuillet C, Quesneville H. Considering transposable element diversification in de novo annotation approaches. PLoS One. 2011 Jan 31;6(1):e16526. doi: 10.1371/journal.pone.0016526. 
* Khaja R, MacDonald JR, Zhang J, Scherer SW. Methods for identifying and mapping recent segmental and gene duplications in eukaryotic genomes. Methods Mol Biol. 2006;338:9-20.
* Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015 <http://www.repeatmasker.org>. 
