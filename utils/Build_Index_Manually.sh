#!/usr/bin/sh
# this script were made for build index in your server as input reference of circPipe. 
# To note , the index building step will take hours and a large amount of resources of your server. You can run the commmand below and have a sleep.  Then start your analysis trip with circPipe 



if [[ $# -eq 0 ]] ; then
    echo 'Usage:\t\tsh Build_index_Manually.sh <genome.fa> <genome.gtf> <basename>'
    echo '\t\t<genome.fa> is your genome sequence in fasta format'
    echo '\t\t<genome.gtf> is the GTF file of your genome sequence'
    echo '\t\t<basename> is your index name. e.g. hg19'
    exit 0
fi

GENOME_FASTA=$1
GENOME_GTF=$2
REFNAME=$3


#### Build STAR INDEX 
#### Parameters : --runThreadN means the cpu numbers , --sjdbGTFfile is the GTF Formatted annotation file , --genomeDir is the path you create for index , --genomeFastaFiles is the FASTA Formatted genome file .
#### The example commandline shows below :
STAR --runMode genomeGenerate --runThreadN 8 --sjdbGTFfile ${GENOME_GTF} --genomeDir ${REFNAME}_starindex/ --genomeFastaFiles ${GENOME_FASTA} --sjdbOverhang 149

#### Build BOWTIE2 INDEX
#### Parameters : ${GENOME_FASTA} is the FASTA Formatted genome file , 'genome' is the prefix name of the index.
bowtie2-build -f ${GENOME_FASTA} ${REFNAME}

#### Build BOWTIE INDEX
#### Parameters : ${GENOME_FASTA} is the FASTA Formatted genome file , 'genome' is the prefix name of the index.
bowtie-build ${GENOME_FASTA} ${REFNAME}

#### Build BWA INDEX
#### Parameters : index '${GENOME_FASTA}' is the FASTA Formatted genome file , -p 'genome' is the prefix name of the index.

bwa index ${GENOME_FASTA} -p ${REFNAME}

#### Build SEGEMEHL INDEX
#### Parameters : -d '${GENOME_FASTA}' is the FASTA Formatted genome file , -x 'genome.idx' is the whole name of the index.
segemehl.x -d ${GENOME_FASTA} -x ${REFNAME}.idx