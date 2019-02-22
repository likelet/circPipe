#!/bin/sh

merge_bed=$1 #merge bed6 file of all piplines results 最终所有样本所有软件的bed6文件
thread=$2 #thread for index buliding # 运行线程
genome=$3 #genome fasta file path

tmp_candidate_hisat_index="candidate_circRNA_index"

# sort bed (in some result bed file , start > end ) and length filtering( >= 100nt)
 awk -F  "\t" '{OFS="\t"}{if ($3 > $2) {name=($1"_"$2"_"$3"_"$6);print $1,$2,$3,name,$5,$6} else {name=($1"_"$3"_"$2"_"$6);print $1,$3,$2,name,$5,$6} }' $merge_bed | awk '$3 - $2 >= 100 ' >  tmp_candidate_circRNA.bed

# bed to gff3 for htseq-count; sites around junction sites(+/-3bp)
awk  '{OFS="\t"}{split($4,a,"_");len=$3-$2; print $4"("a[4]")",".","exon",len-3,len+3,".","+",".","gene_id="$4 }' tmp_candidate_circRNA.bed > tmp_candidate_circRNA.gff3

# bed to fasta
bedtools getfasta -fi $genome -s -bed tmp_candidate_circRNA.bed -name > tmp_candidate.circular.fa

# candidate circRNA sequnces (doulbed).
awk 'NR%2==1{print $0}NR%2==0{print $1$1}'  tmp_candidate.circular.fa > tmp_candidate.circular_doulbed.fa

#build index for candidate circRNA sequnce
mkdir $tmp_candidate_hisat_index

hisat2-build -p $thread tmp_candidate.circular_doulbed.fa $tmp_candidate_hisat_index/candidate_circRNA_doulbed 


