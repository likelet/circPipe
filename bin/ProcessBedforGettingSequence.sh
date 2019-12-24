#!/bin/sh
merge_bed=$1
sortBed=$2 # for bed with less than 400 bp, use all sequence 
startBed=$3
endBed=$4
overhangSize=1000
# sort bed (in some result bed file , start > end ) and length filtering( >= 100nt && <=100kb)
# circRNA > 400 bp
awk -F  "\t" '{OFS="\t"}{if ($3 > $2) {name=($1"_"$2"_"$3"_"$6);print $1,$2,$3,name,$5,$6} else {name=($1"_"$3"_"$2"_"$6);print $1,$3,$2,name,$5,$6} }' $merge_bed | awk -v overhangSize=$overhangSize '$3 - $2 > overhangSize*2 && $3 - $2 <=100000 ' >  tmp_candidate_circRNA.bed
# circRNA <= 400 bp
awk -F  "\t" '{OFS="\t"}{if ($3 > $2) {name=($1"_"$2"_"$3"_"$6);print $1,$2,$3,name,$5,$6} else {name=($1"_"$3"_"$2"_"$6);print $1,$3,$2,name,$5,$6} }' $merge_bed | awk -v overhangSize=$overhangSize '$3 - $2 >= 100 && $3 - $2 <= overhangSize*2 ' >  $sortBed
# star and end site
awk -v overhangSize=$overhangSize '{OFS="\t"}{if ($6=="+"){print $1,$2,$2+overhangSize,$4,"start",$6} else print $1,$3-overhangSize,$3,$4,"start",$6}' tmp_candidate_circRNA.bed > $startBed
awk -v overhangSize=$overhangSize '{OFS="\t"}{if ($6=="+"){print $1,$3-overhangSize,$3,$4,"end",$6} else {print $1,$2,$2+overhangSize,$4,"end",$6} }' tmp_candidate_circRNA.bed > $endBed