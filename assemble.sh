#!/bin/bash
file=rawdata/$1
strain=$1
adapt=~/compbio/adapters/truseq_adapters.fasta

scythe -a $adapt $file.1.fastq -o $file.scythed.1.fastq -q sanger
scythe -a $adapt $file.2.fastq -o $file.scythed.2.fastq -q sanger

sickle pe -f $file.scythed.1.fastq -r $file.scythed.2.fastq -o $file.sickled.1.fastq -p $file.sickled.2.fastq -s $file.singles.fastq -t sanger

pear -f $file.sickled.1.fastq -r $file.sickled.2.fastq -o $file -j 8

cat $file.singles.fastq $file.assembled.fastq > $file.all.singles.fastq

spades.py -m 16 -s $file.all.singles.fastq -1 $file.unassembled.forward.fastq -2 $file.unassembled.reverse.fastq --careful --cov-cutoff auto -t 8 -o $file

## use GetGenomeStats.py to identify appropriate cutoffs for contig length and coverage

## prokka --genus Pseudomonas --strain $strain --locustag $strain --prefix $strain --cpus 8 --outdir $file.prokka $file.filteredcontigs.fna
