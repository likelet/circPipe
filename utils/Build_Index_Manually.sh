#!/usr/bin/sh
# this script were made for build index in your server as input reference of circPipe. 
# To note , the index building step will take hours and a large amount of resources of your server. You can run the commmand below and have a sleep.(Except building the KNIFE index)  Then start your analysis trip with circPipe



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
EXON_GTF=$4
path_to_output_directory=$5
path_to_circularRNApipeline=$6
window_size=$7
transcriptome=$8
REFNAME2=$9
ribosomal=$10
REFNAME3=$11


#### Build STAR INDEX 
#### Parameters : --runThreadN means the cpu numbers , --sjdbGTFfile is the GTF Formatted annotation file , --genomeDir is the path you create for index , --genomeFastaFiles is the FASTA Formatted genome file .
#### The example commandline shows below :
STAR --runMode genomeGenerate --runThreadN 8 --sjdbGTFfile ${GENOME_GTF} --genomeDir ${REFNAME}_starindex/ --genomeFastaFiles ${GENOME_FASTA} --sjdbOverhang 149

#### Build BOWTIE2 INDEX
#### Parameters : ${GENOME_FASTA} is the FASTA Formatted genome file , ${REFNAME} is the prefix name of the index.
bowtie2-build -f ${GENOME_FASTA} ${REFNAME}

#### Build BOWTIE INDEX
#### Parameters : ${GENOME_FASTA} is the FASTA Formatted genome file , ${REFNAME} is the prefix name of the index.
bowtie-build ${GENOME_FASTA} ${REFNAME}

#### Build BWA INDEX
#### Parameters : index '${GENOME_FASTA}' is the FASTA Formatted genome file , -p ${REFNAME} is the prefix name of the index.

bwa index ${GENOME_FASTA} -p ${REFNAME}

#### Build SEGEMEHL INDEX
#### Parameters : -d '${GENOME_FASTA}' is the FASTA Formatted genome file , -x '${REFNAME}.idx' is the whole name of the index.
segemehl.x -d ${GENOME_FASTA} -x ${REFNAME}.idx




####Build KNIFE INDEX
####Please follow the proceduces below manually , or you can take a look at https://github.com/lindaszabo/KNIFE/tree/master/createJunctionIndex
###########################
# Package Requirements
###########################

#python (tested with 2.7.5) and the biopython and bcbio-gff libraries

#Bowtie2.2.2

############################
# Before beginning
############################

#Select 2 string identifiers that will be used for this species index.
#These strings will be used in the subsequent steps and any time you see INDEX_MODE_ID or
#INDEX_FILE_ID in the instructions make sure you use these exact strings.

#    INDEX_MODE_ID: a string that will be included in the MODE parameter when calling findCircularRNA.sh or completeRun.sh to specify alignments to this index
#    INDEX_FILE_ID: a string that will be used in all file names for fasta and bowtie indices for this species

#These can be the same string, but generally it is good to have the file name identifier (#2) be more specific to include information like what genome build and
#annotation are used, while the string for the MODE parameter (#1) is more simple (i.e. mouse, fly, etc).

#############################################################################################
# STEP 1: Gather required data files and create genome, transcriptome, and ribosomal indices
#############################################################################################

#You will need to collect/create the files listed below in order to create an index for a new species
#and run the circular RNA analysis pipeline.

#1. fasta file containing the sequences of the genome
#2. fasta file containing the sequences of the transcriptome
#3. fasta file containing the ribosomal RNA sequences (can be split into subunits or just a single sequence)
#4. a gtf file for the species. Junctions are reported based on the "gene_name" field. This is present in most gtf formats,
#   but if the gtf for your species labels the field containing the gene name something else, you will need to update the
#   gtf to change the label the field to "gene_name" before creating junction indices.
#5. either find or create a bowtie2 genome index that was created with the same genome build you are using for 1-4
#6. either find or create a bowtie genome index. This is bowtie 1, not bowtie2, and is used in the denovo pipeline
#7. either find or create a bowtie2 transcriptome index that was created with the same gtf and genome build you are using

bowtie2-build -f ${transcriptome} ${REFNAME2}

#8. create a bowtie2 ribosomal RNA index using the ribosomal fasta file from #3

bowtie2-build -f ${ribosomal} ${REFNAME3}


#1-6 usually found by downloading igenomes for the species

#7 you can get by downloading a fasta file of known genes from UCSC for the species and using this to create a bowtie2 index.
#Note that when you download the fasta file from UCSC, the gene headers have spaces in them. These need to be removed before
#using the file to create an index, because bowtie will ignore everything after the first space in the header which usually means
#you lose important information. For a fasta file named transcriptome.fa, this command will do the trick:

#more transcriptome.fa | sed -e 's/\s*//g' > transcriptome_fullheaders.fa

#8 takes some searching in the UCSC genome browser or other sources to figure out coordinates of ribosomal subunits.
#If you don't care about subunits, the full ribosomal repeat sequence is usually included in a file downloaded from igenomes
#and you can just use this as your ribosomal RNA.

#######################################################
# STEP 2: Create linear and scrambled junction indices
#######################################################

#There are 2 bash commands you need to call for this. The 1st must be complete before you can run the 2nd. After both are complete, you will find
#the Bowtie2 index files and fasta files for both the linear and scrambled junction indices can be found under circularRNApipeline/index
#and are ready to use.

#COMMAND 1: creates the exon database. This is independent of window size, so if you decide to change the window size you do not need to repeat this step.

python makeExonDB.py -f ${GENOME_FASTA} -a ${EXON_GTF} -o ${path_to_output_directory}

#parameters

#-f: ${GENOME_FASTA} file for the species

#-a: ${EXON_GTF} file for the species with exon annotations to be used for the index

#-o: directory to output the database files

#COMMAND 2: creates the junction indices. This can only be run after COMMAND 1 is complete.

sh createJunctionIndex.sh ${path_to_circularRNApipeline} ${path_to_output_directory} ${REFNAME} ${window_size}

#parameters

#path_to_circularRNApipeline: path to the pipeline directory, index files created will be placed in the index subdirectory

#path_to_output_directory: the same outputDirectory used in COMMAND 1 to create the exonDB

#REFNAME: string selected above that will be included in all of the fasta file names and bowtie index file names

#window_size: (default 1,000,000). size of sliding window to include exon pairs in junction database.

##############################################################################
# STEP 3: Place other index files in correct location to be used by pipeline
##############################################################################

#Rename the genome, transcriptome, and ribosomal files as specified below and place all in circularRNApipeline/index.
#Note that your bowtie2 index file extensions may be .bt2, that is fine, use the existing extension.

#    Bowtie2 genome index: rename all 6 index files to REFNAME_genome.*.bt2l
#    Bowtie2 transcriptome index: rename all 6 index files to REFNAME_transcriptome.*.bt2l
#    Bowtie2 ribosomal index: rename all 6 index files to REFNAME_ribosomal.*.bt2l
#    Fasta files: name the fasta files used to generate these indices respectively to REFNAME_genome.fa, REFNAME_transcriptome.fa, REFNAME_ribosomal.fa

#Rename these files and place in circularRNApipeline/denovo_scripts/index

#    Bowtie1 genome index: rename all 6 index files to REFNAME_genome.*.ebwt

#The gtf file needs to be renamed to REFNAME_genes.gtf and placed in circularRNApipeline/denovo_scripts

#######################################################
# STEP 4: Edit scripts to enable use of the new index
#######################################################

#Annotation-dependent script: update circularRNApipeline/findCircularRNA.sh to add a bt_prefix condition for your species.
#bt_prefix will be set to REFNAME and the mode string to search for will be the REFNAME

#Denovo script: update circularRNApipeline/denovo_scripts/process_directory_unaligned.pl to add a case to point $reference and $gtf variables
#to your files when $mode contains REFNAME