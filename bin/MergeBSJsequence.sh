shortFa=$1
startFa=$2
endFa=$3
outFa=$4


### BSJ flanking region sequence (+/- 200bp).
awk 'NR==FNR{if(NR%2==1) {getline $2;a[$1]=$2 }}NR>FNR {if(NR%2==1){getline $2;{print $1"\n"$2a[$1]}}}'  $startFa $endFa > tmp_candidate.circular_BSJ_flank1.fa

awk 'NR%2==1{print $0}NR%2==0{print $1$1}' $shortFa > tmp_candidate.circular_BSJ_flank2.fa

cat tmp_candidate.circular_BSJ_flank1.fa tmp_candidate.circular_BSJ_flank2.fa > $outFa