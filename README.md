# SLHSD
SLHSD is a Hybrid scaffolding method based on short and long reads
=========


Scaffolder: SLHSD
=================

1) Introduction
```
    SLHSD is an scaffolder which aims to determine the orientations and orders of contigs. 
    The contigs can be produced by any assembler.
    The input data of SLHSD is the short reads (fasta format),long reads (fasta format) and the contigs (fasta format). 
```
2) Before installing and running
```
    Please install BWA from https://github.com/lh3/bwa.
	Please install Samtools from https://sourceforge.net/projects/samtools/files/samtools/.
	Please build and install Bamtools from https://github.com/pezmaster31/bamtools.
```
3) Installing.
```
    SLHSD should run on Linux operating sysetm with gcc. We test SLHSD using gcc9.4.0 on Ubuntu.
    Create a main directory. Copy all source code to this directory.
	cd SLHSD
	export BAMTOOLS_HOME_INCLUDE=/path_bamtools_include_api_shared/
	export BAMTOOLS_HOME_LIB=/path_bamtools_lib_libbamtools.a/
	make all
```
4) Running.
```
	step 1: bwa index contigs.fasta
	step 2: bwa mem -a contigs.fasta shortreads1.fasta > sr-align1.sam
	step 3: samtools view -Sb sr-align1.sam > sr-align1.bam
	step 4: bwa mem -a contigs.fasta shortreads2.fasta > sr-align2.sam
	step 5: samtools view -Sb sr-align2.sam > sr-align2.bam
	step 6: bwa mem -t8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y contigs.fasta longreads.fasta > lr-align.sam
	step 7: samtools view -Sb lr-align.sam > lr-align.bam
	step 8: SLHSD -c contigs.fasta -l sr-align1.bam -r sr-align2.bam -b lr-align.bam -p short-read-type -d insertSize
	
	-c <contigs.fa>: 
	    The file includes contigs produced by one assembler.
	-l <sr-align1.bam>:
	    The aligning result between the contigs and the short left reads.
	-r <sr-align2.bam>:
	    The aligning result between the contigs and the short right reads.
	-b <lr-align3.sam>:
	    The aligning result between the contigs and the long reads.
	-p <readType>: 
	    The type of short reads: paired-end(p=1) or mate-pair(p=2).
	-d <insertSize>: 
	    The insert size of read library.
	-a <readLength>: 
		The length of read.
			
```
5) Output.
```
    The output file "scaffold_set.fa" is the scaffolding result. 
```


