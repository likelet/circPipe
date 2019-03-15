#!/usr/bin/sh

#### STAR INDEX
#### Parameters : --runThreadN means the cpu numbers , --sjdbGTFfile is the GTF Formatted annotation file , --genomeDir is the path you create for index , --genomeFastaFiles is the FASTA Formatted genome file .
#### The example commandline shows below :
STAR --runMode genomeGenerate --runThreadN 8 --sjdbGTFfile gencode.v25.annotation.gtf --genomeDir starindex/ --genomeFastaFiles genome_hg38.fa --sjdbOverhang 149

#### BOWTIE2 INDEX
#### Parameters : genome_hg38.fa is the FASTA Formatted genome file , 'genome' is the prefix name of the index.
#### The example commandline shows below :
bowtie2-build -f genome_hg38.fa genome

#### BOWTIE INDEX
#### Parameters : genome_hg38.fa is the FASTA Formatted genome file , 'genome' is the prefix name of the index.
#### The example commandline shows below :
bowtie-build genome_hg38.fa genome

#### BWA INDEX
#### Parameters : index 'genome_hg38.fa' is the FASTA Formatted genome file , -p 'genome' is the prefix name of the index.
#### The example commandline shows below :
bwa index genome_hg38.fa -p genome

#### SEGEMEHL INDEX
#### Parameters : -d 'genome_hg38.fa' is the FASTA Formatted genome file , -x 'genome.idx' is the whole name of the index.
#### The example commandline shows below :
segemehl.x -d genome_hg38.fa -x genome.idx