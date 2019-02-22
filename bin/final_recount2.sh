#!/bin/sh
###$4=fastq $5=fastq2
id=$1
seq_type=$2 ###single_end , paired_end
thread=$3

gff="tmp_candidate_circRNA.gff3"
index="candidate_circRNA_index/candidate_circRNA_doulbed"


if [[ "${seq_type}" = "single_end" ]] ; then
	#statements
	hisat2 -p $thread -t -k 1 -x $index -U $4 | samtools view -bS - > $id_circRNA.bam 
elif [[ "${seq_type}" = "paired_end" ]]; then
	fastq1=$4
	fastq2=$5
	hisat2 -p $thread -t -k 1 --no-mixed -x $index -1 $fastq1 -2 $fastq2 | samtools view -bS - > ${id}_circRNA.bam 
fi

#-r name for paired end data
htseq-count -f bam  -s no  $id_circRNA.bam  $gff  > ${id}_circRNA_requantity.count 
