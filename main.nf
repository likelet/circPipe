#!/usr/bin/env nextflow

/*
========================================================================================
                            cirPipe
========================================================================================
 * cirPipe was implemented by Dr. Qi Zhao and Qijin Wei from Sun Yat-sen University Cancer Center.
 * Homepage / Documentation
  https://github.com/likelet/cirpipe

 */

/*
========================================================================================
                         nf-core/cirpipe
========================================================================================
 nf-core/cirpipe Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/cirpipe
----------------------------------------------------------------------------------------
*/


/*
 * to be added
 *
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>: design and implement the pipeline.
 * Wei qijin
 */



//pre-defined functions for render command
//=======================================================================================
ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";
def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }




def helpMessage() {
    log.info"""
    =========================================
     nf-core/cirpipe v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/cirpipe --reads '*_R{1,2}.fastq.gz' -profile standard,docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

/*
// Configurable variables
params.name = false
//params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
//params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

output_docs = file("$baseDir/docs/output.md")


// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

 *
 * Create a channel for input read files

if(params.readPaths){
    if(params.singleEnd){
        Channel
                .from(params.readPaths)
                .map { row -> [ row[0], [file(row[1][0])]] }
                .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
                .into { read_files_fastqc; read_files_trimming }
    } else {
        Channel
                .from(params.readPaths)
                .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
                .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
                .into { read_files_fastqc; read_files_trimming }
    }
} else {
    Channel
            .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .into { read_files_fastqc; read_files_trimming }
}

*/

/*
// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/cirpipe v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/cirpipe'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile

/*
if(workflow.profile == 'awsbatch'){
    summary['AWS Region'] = params.awsregion
    summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="



 *
 * Checking the input files
 * Build up the output files
 * Adding input files error exceptions Here
 */

outdir = file(params.outdir) //the output directory
if( !outdir.exists() ) exit 1, print_red("Missing output directory: ${output}")

def fastpoutdir = new File( "${params.outdir}/pipeline_fastp/" )
if( !fastpoutdir.exists() ) {
    fastpoutdir.mkdirs()
}

def toolsoutdir = new File( "${params.outdir}/pipeline_tools/" )
if( !toolsoutdir.exists() ) {
    toolsoutdir.mkdirs()
}

def staroutdir = new File( "${params.outdir}/pipeline_tools/pipeline_star/" )
if( !staroutdir.exists() ) {
    staroutdir.mkdirs()
}

def bwaoutdir = new File( "${params.outdir}/pipeline_tools/pipeline_bwa/" )
if( !bwaoutdir.exists() ) {
    bwaoutdir.mkdirs()
}

def bowtie2outdir = new File( "${params.outdir}/pipeline_tools/pipeline_bowtie2/" )
if( !bowtie2outdir.exists() ) {
    bowtie2outdir.mkdirs()
}

def circexplorer2outdir = new File( "${params.outdir}/pipeline_tools/pipeline_circexplorer2/" )
if( !circexplorer2outdir.exists() ) {
    circexplorer2outdir.mkdirs()
}

def cirioutdir = new File( "${params.outdir}/pipeline_tools/pipeline_ciri/" )
if( !cirioutdir.exists() ) {
    cirioutdir.mkdirs()
}

def find_circoutdir = new File( "${params.outdir}/pipeline_tools/pipeline_find_circ/" )
if( !find_circoutdir.exists() ) {
    find_circoutdir.mkdirs()
}

def mapspliceoutdir = new File( "${params.outdir}/pipeline_tools/pipeline_mapsplice/" )
if( !mapspliceoutdir.exists() ) {
    mapspliceoutdir.mkdirs()
}

def segemehloutdir = new File( "${params.outdir}/pipeline_tools/pipeline_segemehl/" )
if( !segemehloutdir.exists() ) {
    segemehloutdir.mkdirs()
}

def mergeoutdir = new File( "${params.outdir}/pipeline_tools/pipeline_merge/" )
if( !mergeoutdir.exists() ) {
    mergeoutdir.mkdirs()
}

def DEoutdir = new File( "${params.outdir}/DE/" )
if( !DEoutdir.exists() ) {
    DEoutdir.mkdirs()
}

def CORoutdir = new File( "${params.outdir}/COR/" )
if( !CORoutdir.exists() ) {
    CORoutdir.mkdirs()
}

def REPORToutdir = new File( "${params.outdir}/REPORT/" )
if( !REPORToutdir.exists() ) {
    REPORToutdir.mkdirs()
}

otherTools = file(params.otherTools) //the other tools directory
if( !otherTools.exists() ) exit 1, print_red("Missing other tools directory: ${otherTools}")


starindex = file(params.starindex) //the index directory
if( !starindex.exists() ) exit 1, print_red("Missing star index directory: ${starindex}")

segindex = file(params.segindex) //the segemehl index file
if( !segindex.exists() ) exit 1, print_red("Missing Segemehl index file: ${segindex}")


refdir = file(params.refdir) //the reference genome directory
if( !refdir.exists() ) exit 1, print_red("Missing Reference Genome Directory: ${refdir}")

annotationfile = file(params.annotationfile) //the annotationfile
if( !annotationfile.exists() ) exit 1, print_red("Missing annotation file: ${annotationfile}")

genomefile = file(params.genomefile) //the genomefile
if( !genomefile.exists() ) exit 1, print_red("Missing genome file: ${genomefile}")

gtffile = file(params.gtffile) //the annotationfile
if( !gtffile.exists() ) exit 1, print_red("Missing annotation file: ${gtffile}")



condadir = file(params.condadir) //the python3 environment
if( !condadir.exists() ) exit 1, print_red("Missing python3 environment: ${condadir}")

conda2dir = file(params.conda2dir) //the python2 environment
if( !conda2dir.exists() ) exit 1, print_red("Missing python2 environment: ${conda2dir}")


mapsdir = file(params.mapsdir) //the mapsplice directory
if( !mapsdir.exists() ) exit 1, print_red("Missing Mapsplice Directory: ${mapsdir}")

segdir = file(params.mapsdir) //the segemehl directory
if( !segdir.exists() ) exit 1, print_red("Missing Segemehl Directory: ${segdir}")

ciridir = file(params.ciridir)
if( !genomefile.exists() ) exit 1, print_red("Missing CIRI Directory: ${ciridir}")

find_circdir = file(params.find_circdir)
if( !find_circdir.exists() ) exit 1, print_red("Missing find_circ Directory: ${find_circdir}")




//showing the process and files
log.info print_purple("""\
         c i r P i p e   P I P E L I N E
         =============================


         Reads types :
         singleEnd : ${params.singleEnd}
         
         Tools selected :
         star : ${params.star}
         bwa : ${params.bwa}
         bowtie2 : ${params.bowtie2}
         mapsplice : ${params.mapsplice}
         segemehl : ${params.segemehl}
         selectAll : ${params.selectAll}

         Input files selected :
         reads : ${params.reads}
         annotation file : ${params.annotationfile}
         genome file : ${params.genomefile}
         gtf file : ${params.gtffile}

         Output files directory :
         output directory : ${params.outdir}


         Start running...


         """)
        .stripIndent()



ava_cpu = Runtime.getRuntime().availableProcessors()
// set individual cpu for fork run
if ( params.cpu != null && ava_cpu > params.cpu ) {
    idv_cpu = params.cpu
} else if ( params.cpu != null && ava_cpu < params.cpu ) {
    idv_cpu = ava_cpu
    print print_yellow("Exceeding the max available processors, \n use default parameter to run pipe. ")
}
int fork_number = ava_cpu / idv_cpu
if (fork_number < 1) {
    fork_number = 1
}

/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
if(params.singleEnd){
    Channel
            .fromFilePairs( params.readssingleEnd, size: 1 )
            .ifEmpty { error "Cannot find any reads matching: ${params.readssingleEnd}" }
            .set { read_pairs_fastp }
}else{
    Channel
            .fromFilePairs( params.reads )
            .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
            .set { read_pairs_fastp }
}



//run the fastp
process run_fastp{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_fastp", mode: 'copy', pattern: "*_fastpreport.html", overwrite: true

    maxForks fork_number

    input:
    set pair_id, file(query_file) from read_pairs_fastp

    output:
    set pair_id, file ('fastp_*') into fastpfiles_star
    set pair_id, file ('fastp_*') into fastpfiles_bwa
    set pair_id, file ('unzip_fastp_*') into fastpfiles_mapsplice
    set pair_id, file ('unzip_fastp_*') into fastpfiles_segemehl
    set pair_id, file ('fastp_*') into fastpfiles_bowtie2
    file ('*.html') into fastp_for_waiting
    file ('*_fastp.json') into fastp_for_multiqc

    script:
    if(params.singleEnd){
        """
        fastp \
        -i ${query_file} \
        -o unzip_fastp_${pair_id}.fq

        fastp \
        -i ${query_file} \
        -o fastp_${pair_id}.fq.gz \
        -h ${pair_id}_fastpreport.html \
        -j ${pair_id}_fastp.json
        """
    }else{
        """
        fastp \
        -i ${query_file[0]} \
        -I ${query_file[1]} \
        -o unzip_fastp_${pair_id}_1.fq \
        -O unzip_fastp_${pair_id}_2.fq

        fastp \
        -i ${query_file[0]} \
        -I ${query_file[1]} \
        -o fastp_${pair_id}_1.fq.gz \
        -O fastp_${pair_id}_2.fq.gz \
        -h ${pair_id}_fastpreport.html \
        -j ${pair_id}_fastp.json
        """
    }


}

fastp_for_waiting = fastp_for_waiting.first() //wait for finish this process first

//run the multiqc
process run_multiqc{
    publishDir "${params.outdir}/pipeline_fastp", mode: 'copy', pattern: "*.html", overwrite: true

    maxForks fork_number

    input:
    file (query_file) from fastp_for_multiqc.collect()

    output:
    file ('*.html') into multiqc_results

    script:
    """
    multiqc .
    """
}


//the first tool : star - circexplorer2
//run the star
process run_star{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_star", mode: 'link', overwrite: true

    maxForks fork_number

    input:
    set pair_id, file(query_file) from fastpfiles_star
    file starindex

    output:
    set pair_id, file ('*.junction') into starfiles

    when:
    params.circexplorer2 || params.selectAll

    shell:
    star_threads = idv_cpu - 1
    if(params.singleEnd){
        """
        STAR \
        --runThreadN ${task.cpus} \
        --chimSegmentMin 10 \
        --genomeDir ${starindex} \
        --readFilesCommand zcat \
        --readFilesIn ${query_file} \
        --outFileNamePrefix star_${pair_id}_
        """
    }else{
        """
        STAR \
        --runThreadN ${task.cpus} \
        --chimSegmentMin 10 \
        --genomeDir ${starindex} \
        --readFilesCommand zcat \
        --readFilesIn ${query_file[0]} ${query_file[1]} \
        --outFileNamePrefix star_${pair_id}_
        """
    }

}

//run the circexplorer2
process run_circexplorer2{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_circexplorer2", mode: 'copy', overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from starfiles
    file annotationfile
    file genomefile

    output:
    set pair_id, file ('*known.txt') into circexplorer2files

    when:
    params.circexplorer2 || params.selectAll

    script:
    """
        CIRCexplorer2 \
        parse -t STAR ${query_file} \
        > CIRCexplorer2_parse_${pair_id}.log

        CIRCexplorer2 \
        annotate -r ${annotationfile} \
        -g ${genomefile} \
        -b back_spliced_junction.bed \
        -o CIRCexplorer2_${pair_id}_circularRNA_known.txt \
        > CIRCexplorer2_annotate_${pair_id}.log
        """
}

//produce the bed6 file
process run_modify_circexplorer2{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_circexplorer2", mode: 'copy', pattern: "*candidates.bed", overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from circexplorer2files
    file otherTools

    output:
    file ('*candidates.bed') into modify_circexplorer2
    val (pair_id) into modify_circexplorer2_id

    when:
    params.circexplorer2 || params.selectAll

    shell :
    '''
        grep circ !{query_file} \
        | grep -v chrM \
	    | awk '{print $1 "\t" $2 "\t" $3 "\t" "circexplorer2" "\t" $13 "\t" $6}' \
        > !{pair_id}_modify_circexplorer2.temp.bed
        
        python !{otherTools}/quchong.py !{pair_id}_modify_circexplorer2.temp.bed !{pair_id}_modify_circexplorer2.candidates.bed
        '''
}

//produce the matrix
process matrix_circexplorer2{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_merge", mode: 'copy', pattern: "*.matrix", overwrite: true

    maxForks fork_number

    input:
    file (query_file) from modify_circexplorer2.collect()
    val (pair_id) from modify_circexplorer2_id.collect()
    file otherTools

    output:
    file ('circexplorer2.txt') into merge_circexplorer2
    file ('*.matrix') into output_circexplorer2
    file ('name_circexplorer2.txt') into name_circexplorer2

    when:
    params.circexplorer2 || params.selectAll

    shell :
    '''
		for file in !{query_file}
		do
			cat $file >> concatenate.bed
		done
		
		python !{otherTools}/hebinglist.py concatenate.bed merge_concatenate.bed
		sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
		cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
		cat mergeconcatenate.bed > circexplorer2.txt
		
		for file in !{query_file}
		do
			python !{otherTools}/quchongsamples.py mergeconcatenate.bed $file counts.txt
			paste -d"\t" id.txt counts.txt > temp.txt
			cat temp.txt > id.txt
		done
		
		echo -e "id\\c" > merge_header.txt
        for sampleid in !{pair_id}
        do
            echo -e "\\t$sampleid\\c" >> merge_header.txt    
        done 
        
        sed -i 's/\\[//g' merge_header.txt
        sed -i 's/\\,//g' merge_header.txt
        sed -i 's/\\]//g' merge_header.txt
        echo -e "\\n\\c" >> merge_header.txt
         
        cat merge_header.txt id.txt > circexplorer2_merge.matrix
        echo -e "circexplorer2" > name_circexplorer2.txt
        '''
}


//the second tool : bwa - ciri
//run the bwa
process run_bwa{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_bwa", mode: 'link', overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from fastpfiles_bwa

    output:
    set pair_id, file ('*.sam') into bwafiles

    when:
    params.ciri || params.selectAll

    shell:
    bwa_threads = idv_cpu - 1
    if(params.singleEnd){
        """
        bwa \
        mem -t ${task.cpus} \
        -k 15 \
        -T 19  -M -R \
        "@RG\\tID:fastp_${pair_id}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${pair_id}" \
        /home/wqj/test/bwaindex/genome \
        ${query_file} \
        > bwa_${pair_id}.mem.sam
        """
    }else{
        """
        bwa \
        mem -t ${task.cpus} \
        -T 19 -M -R \
        "@RG\\tID:fastp_${pair_id}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${pair_id}" \
        /home/wqj/test/bwaindex/genome \
        ${query_file[0]} ${query_file[1]} \
        > bwa_${pair_id}.mem.sam
        """
    }

}

//run the ciri
process run_ciri{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_ciri", mode: 'copy', overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from bwafiles
    file gtffile
    file genomefile
    file ciridir

    output:
    set pair_id, file ('*.txt') into cirifiles

    when:
    params.ciri || params.selectAll

    script:
    """
        perl ${ciridir}/CIRI2.pl \
        -T 10 \
        -F ${genomefile} \
        -A ${gtffile} \
        -G CIRI_${pair_id}.log \
        -I ${query_file} \
        -O CIRI_${pair_id}.txt \
        > CIRI_${pair_id}_detail.log
        """
}

//produce the bed6 file
process run_modify_ciri{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_ciri", mode: 'copy', pattern: "*candidates.bed", overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from cirifiles
    file otherTools

    output:
    file ('*candidates.bed') into modify_ciri_file
    val (pair_id) into modify_ciri_id

    when:
    params.ciri || params.selectAll

    shell :
    '''
        cat !{query_file} \
	    | sed -e '1d' \
        | grep -v chrM \
        | awk '{print $2 "\t" $3 "\t" $4 "\t" "ciri" "\t" $5 "\t" $11}' \
        > !{pair_id}_modify_ciri.temp.bed
        
        python !{otherTools}/quchong.py !{pair_id}_modify_ciri.temp.bed !{pair_id}_modify_ciri.candidates.bed
        '''
}

//produce the matrix
process matrix_ciri{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_merge", mode: 'copy', pattern: "*.matrix", overwrite: true

    maxForks fork_number

    input:
    file (query_file) from modify_ciri_file.collect()
    val (pair_id) from modify_ciri_id.collect()
    file otherTools

    output:
    file ('ciri.txt') into merge_ciri
    file ('*.matrix') into output_ciri
    file ('name_ciri.txt') into name_ciri

    when:
    params.ciri || params.selectAll

    shell :
    '''
		for file in !{query_file}
		do
			cat $file >> concatenate.bed
		done
		
		python !{otherTools}/hebinglist.py concatenate.bed merge_concatenate.bed
		sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
		cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
		cat mergeconcatenate.bed > ciri.txt
		
		for file in !{query_file}
		do
			python !{otherTools}/quchongsamples.py mergeconcatenate.bed $file counts.txt
			paste -d"\t" id.txt counts.txt > temp.txt
			cat temp.txt > id.txt
		done
		
		echo -e "id\\c" > merge_header.txt
        for sampleid in !{pair_id}
        do
            echo -e "\\t$sampleid\\c" >> merge_header.txt    
        done 
        
        sed -i 's/\\[//g' merge_header.txt
        sed -i 's/\\,//g' merge_header.txt
        sed -i 's/\\]//g' merge_header.txt
        echo -e "\\n\\c" >> merge_header.txt
         
        cat merge_header.txt id.txt > ciri_merge.matrix
        echo -e "ciri" > name_ciri.txt
        '''
}


//the third tool : mapsplice
//run the mapsplice
process run_mapsplice{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_mapsplice", mode: 'copy', overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from fastpfiles_mapsplice
    file mapsdir
    file gtffile
    file refdir
    file outdir

    output:
    set pair_id, file('*.log') into mapsplicefiles

    conda params.condadir

    when:
    params.mapsplice || params.selectAll

    shell:
    mapsplice_threads = idv_cpu - 1
    if(params.singleEnd){
        """
        python ${mapsdir}/mapsplice.py \
        -p ${task.cpus} \
        -k 1 \
        --fusion-non-canonical \
        --qual-scale phred33 \
        --non-canonical-double-anchor \
        --min-fusion-distance 200 \
        -x /home/wqj/test/bowtieindex/chrX \
        --gene-gtf ${gtffile} \
        -c ${refdir} \
        -1 ${query_file} \
        -o ${outdir}/pipeline_tools/pipeline_mapsplice/output_mapsplice_${pair_id} 2 \
        > ${pair_id}_mapsplice.log
        """
    }else{
        """
        python ${mapsdir}/mapsplice.py \
        -p ${task.cpus} \
        -k 1 \
        --fusion-non-canonical \
        --non-canonical-double-anchor \
        --min-fusion-distance 200 \
        -x /home/wqj/test/bowtieindex/chrX \
        --gene-gtf ${gtffile} \
        -c ${refdir} \
        -1 ${query_file[0]} \
        -2 ${query_file[1]} \
        -o ${outdir}/pipeline_tools/pipeline_mapsplice/output_mapsplice_${pair_id} 2 \
        > ${pair_id}_mapsplice.log
        """
    }

}

//produce the bed6 file
process run_modify_mapsplice{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_mapsplice", mode: 'copy', pattern: "*candidates.bed", overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from mapsplicefiles
    file outdir
    file otherTools

    output:
    file ('*candidates.bed') into modify_mapsplice
    val (pair_id) into modify_mapsplice_id

    when:
    params.mapsplice || params.selectAll

    shell :
    '''
    cat !{outdir}/pipeline_tools/pipeline_mapsplice/output_mapsplice_!{pair_id}/circular_RNAs.txt \
    | awk '{print $6}' \
    | sed -e 's/.//' \
    > !{pair_id}_mapsplice_temp1.bed

    cat !{outdir}/pipeline_tools/pipeline_mapsplice/output_mapsplice_!{pair_id}/circular_RNAs.txt \
    | awk '{print $1}' \
    | awk -F"~" '{print $2}' \
    > !{pair_id}_mapsplice_temp.bed

    paste !{pair_id}_mapsplice_temp.bed !{pair_id}_mapsplice_temp1.bed !{outdir}/pipeline_tools/pipeline_mapsplice/output_mapsplice_!{pair_id}/circular_RNAs.txt \
    | grep -v chrM \
    | awk '{if($2=="-") print $1 "\t" $4 "\t" $5 "\t" "mapsplice" "\t" $7 "\t" $2 ; else print $1 "\t" $5 "\t" $4 "\t" "mapsplice" "\t" $7 "\t" $2 }' \
    > !{pair_id}_modify_mapsplice.temp.bed
    
    python !{otherTools}/quchong.py !{pair_id}_modify_mapsplice.temp.bed !{pair_id}_modify_mapsplice.candidates.bed
    '''
}

//produce the matrix
process matrix_mapsplice{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_merge", mode: 'copy', pattern: "*.matrix", overwrite: true

    maxForks fork_number

    input:
    file (query_file) from modify_mapsplice.collect()
    val (pair_id) from modify_mapsplice_id.collect()
    file otherTools

    output:
    file ('mapsplice.txt') into merge_mapsplice
    file ('*.matrix') into output_mapsplice
    file ('name_mapsplice.txt') into name_mapsplice

    when:
    params.mapsplice || params.selectAll

    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> concatenate.bed
    done
    
    python !{otherTools}/hebinglist.py concatenate.bed merge_concatenate.bed
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
    cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
    cat mergeconcatenate.bed > mapsplice.txt
    
    for file in !{query_file}
    do
        python !{otherTools}/quchongsamples.py mergeconcatenate.bed $file counts.txt
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
    done
    
    echo -e "id\\c" > merge_header.txt
    for sampleid in !{pair_id}
    do
        echo -e "\\t$sampleid\\c" >> merge_header.txt    
    done 
    
    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > mapsplice_merge.matrix
    echo -e "mapsplice" > name_mapsplice.txt
    '''
}


//the fourth tool : segemehl
//run the segemehl
process run_segemehl{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_segemehl", mode: 'copy', overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from fastpfiles_segemehl
    file segdir
    file genomefile
    file segindex

    output:
    set pair_id, file ('*splicesites.bed') into segemehlfiles

    conda params.condadir

    when:
    params.segemehl || params.selectAll

    shell:
    segemehl_threads = idv_cpu - 1
    if(params.singleEnd){
        """
        ${segdir}/segemehl.x \
        -d ${genomefile} \
        -i ${segindex} \
        -q ${query_file} \
        -t ${task.cpus} \
        -S \
        | samtools view -bS - \
        | samtools sort -o - deleteme \
        | samtools view -h - \
        > segemehl_${pair_id}_mapped.sam

        ${segdir}/testrealign.x \
        -d ${genomefile} \
        -q segemehl_${pair_id}_mapped.sam \
        -n \
        -U segemehl_${pair_id}_splicesites.bed \
        -T segemehl_${pair_id}_transrealigned.bed
        """
    }else{
        """
        ${segdir}/segemehl.x \
        -d ${genomefile} \
        -i ${segindex} \
        -q ${query_file[0]} \
        -p ${query_file[1]} \
        -t ${task.cpus} \
        -S \
        | samtools view -bS - \
        | samtools sort -o - deleteme \
        | samtools view -h - \
        > segemehl_${pair_id}_mapped.sam

        ${segdir}/testrealign.x \
        -d ${genomefile} \
        -q segemehl_${pair_id}_mapped.sam \
        -n \
        -U segemehl_${pair_id}_splicesites.bed \
        -T segemehl_${pair_id}_transrealigned.bed
        """
    }

}

//produce the bed6 file
process run_modify_segemehl{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_segemehl", mode: 'copy', pattern:"*candidates.bed", overwrite: true

    maxForks fork_number

    input:
    set pair_id , file ( query_file ) from segemehlfiles
    file otherTools

    output:
    file ('*candidates.bed') into modify_segemehl
    val (pair_id) into modify_segemehl_id

    when:
    params.segemehl || params.selectAll

    shell :
    '''
    cat !{query_file} \
    | awk '{print $4}' \
    | awk -F":" '{print $2 "\t" $5 "\t" $6}' \
    > !{pair_id}_segemehl.temp.bed

    paste !{query_file} !{pair_id}_segemehl.temp.bed \
    | grep C \
    | grep P \
    | grep -v chrM \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" "segemehl" "\t" $7 "\t" $6}' \
    > !{pair_id}_modify_segemehl.temp.bed
    
    python !{otherTools}/quchong.py !{pair_id}_modify_segemehl.temp.bed !{pair_id}_modify_segemehl.candidates.bed
    '''
}

//produce the matrix
process matrix_segemehl{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_merge", mode: 'copy', pattern: "*.matrix", overwrite: true

    maxForks fork_number

    input:
    file (query_file) from modify_segemehl.collect()
    val (pair_id) from modify_segemehl_id.collect()
    file otherTools

    output:
    file ('segemehl.txt') into merge_segemehl
    file ('*.matrix') into output_segemehl
    file ('name_segemehl.txt') into name_segemehl

    when:
    params.segemehl || params.selectAll

    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> concatenate.bed
    done
    
    python !{otherTools}/hebinglist.py concatenate.bed merge_concatenate.bed
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
    cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
    cat mergeconcatenate.bed > segemehl.txt
    
    for file in !{query_file}
    do
        python !{otherTools}/quchongsamples.py mergeconcatenate.bed $file counts.txt
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
    done
    
    echo -e "id\\c" > merge_header.txt
    for sampleid in !{pair_id}
    do
        echo -e "\\t$sampleid\\c" >> merge_header.txt    
    done 
    
    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > segemehl_merge.matrix
    echo -e "segemehl" > name_segemehl.txt
    '''
}


//the fifth tool : bowtie2 - find_circ
//run the bowtie2
process run_bowtie2{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_bowtie2", mode: 'link', overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from fastpfiles_bowtie2

    output:
    set pair_id, file ('bowtie2*') into bowtie2files

    when:
    params.find_circ || params.selectAll

    shell:
    bowtie2_threads = idv_cpu - 1
    if(params.singleEnd){
        """
        bowtie2 \
        -p ${task.cpus} \
        --very-sensitive \
        --score-min=C,-15,0 \
        --mm \
        -x /home/wqj/test/bowtie2index/chrX \
        -q \
        -U ${query_file} \
        | samtools view -hbuS - \
        | samtools sort \
        -o bowtie2_output_${pair_id}.bam

        samtools \
        view -hf 4 bowtie2_output_${pair_id}.bam \
        | samtools view -Sb - \
        > bowtie2_unmapped_${pair_id}.bam
        """
    }else{
        """
        bowtie2 \
        -p ${task.cpus} \
        --very-sensitive \
        --score-min=C,-15,0 \
        --mm \
        -x /home/wqj/test/bowtie2index/chrX \
        -q \
        -1 ${query_file[0]} \
        -2 ${query_file[1]} \
        | samtools view -hbuS - \
        | samtools sort \
        -o bowtie2_output_${pair_id}.bam

        samtools \
        view -hf 4 bowtie2_output_${pair_id}.bam \
        | samtools view -Sb - \
        > bowtie2_unmapped_${pair_id}.bam
        """
    }

}

//run the find_circ
process run_find_circ{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_find_circ", mode: 'copy', overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from bowtie2files
    file genomefile
    file find_circdir

    output:
    set pair_id, file ('*splice_sites.bed') into find_circfiles

    conda params.conda2dir

    when:
    params.find_circ || params.selectAll

    shell:
    bowtie2_threads = idv_cpu - 1
    """     
    python ${find_circdir}/unmapped2anchors.py ${query_file} \
    | gzip \
    > find_circ_${pair_id}_anchors.qfa.gz

    bowtie2 \
    -p ${task.cpus} \
    --reorder \
    --mm \
    --score-min=C,-15,0 \
    -q \
    -x /home/wqj/test/bowtie2index/chrX \
    -U find_circ_${pair_id}_anchors.qfa.gz \
    | python ${find_circdir}/find_circ.py \
    -G ${genomefile} \
    -p ${pair_id}_ \
    -s find_circ_${pair_id}_stats.sites.log \
    -n find_circ \
    -R find_circ_${pair_id}_spliced_reads.fa \
    > find_circ_${pair_id}_splice_sites.bed   
    """
}

//produce the bed6 file
process run_modify_find_circ{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_find_circ", mode: 'copy', pattern: "*candidates.bed", overwrite: true

    maxForks fork_number

    input:
    set pair_id, file (query_file) from find_circfiles
    file otherTools

    output:
    file ('*candidates.bed') into modify_find_circfiles
    val (pair_id) into modify_find_circ_id

    when:
    params.bowtie2 || params.selectAll

    shell :
    '''
    grep CIRCULAR !{query_file} \
    | grep -v chrM \
    | grep UNAMBIGUOUS_BP \
    | grep ANCHOR_UNIQUE \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $5 "\t" $6}' \
    | sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 \
    > !{pair_id}_modify_finc_circ.temp.bed
    
    python !{otherTools}/quchong.py !{pair_id}_modify_finc_circ.temp.bed !{pair_id}_modify_finc_circ.candidates.bed
    '''
}

//produce the matrix
process matrix_find_circ{
    tag "$pair_id"
    publishDir "${params.outdir}/pipeline_tools/pipeline_merge", mode: 'copy', pattern: "*.matrix", overwrite: true

    maxForks fork_number

    input:
    file (query_file) from modify_find_circfiles.collect()
    val (pair_id) from modify_find_circ_id.collect()
    file otherTools

    output:
    file ('find_circ.txt') into merge_find_circ
    file ('*.matrix') into output_find_circ
    file ('name_find_circ.txt') into name_find_circ

    when:
    params.find_circ || params.selectAll

    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> concatenate.bed
    done
    
    python !{otherTools}/hebinglist.py concatenate.bed merge_concatenate.bed
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
    cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
    cat mergeconcatenate.bed > find_circ.txt
    
    for file in !{query_file}
    do
        python !{otherTools}/quchongsamples.py mergeconcatenate.bed $file counts.txt
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
    done
    
    echo -e "id\\c" > merge_header.txt
    for sampleid in !{pair_id}
    do
        echo -e "\\t$sampleid\\c" >> merge_header.txt    
    done 
    
    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > find_circ_merge.matrix
    echo -e "find_circ" > name_find_circ.txt
    '''
}


//after running the tools
//calculate the results by different tools
process calculate_tools{
    publishDir "${params.outdir}/pipeline_tools/pipeline_merge", mode: 'copy', pattern: "*.matrix", overwrite: true

    maxForks fork_number

    input:
    file (query_file) from merge_find_circ.concat( merge_circexplorer2, merge_ciri, merge_mapsplice, merge_segemehl ).collect()
    file (name_file) from name_find_circ.concat( name_circexplorer2, name_ciri, name_mapsplice, name_segemehl ).collect()
    file otherTools

    output:
    file ('*.matrix') into tools_merge

    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> concatenate.txt
    done
    
    python !{otherTools}/hebingtoolsid.py concatenate.txt id_unsort.txt
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 id_unsort.txt > id_sort.txt
    
    echo -e "total\\c" > total.txt
    cat id_sort.txt | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id_merge.txt
    cat id_merge.txt total.txt > id_list.txt
    
    for file in !{query_file}
    do
        python !{otherTools}/countnumbers.py id_sort.txt $file counts.txt
        paste -d"\t" id_list.txt counts.txt > temp.txt
        cat temp.txt > id_list.txt
    done
    
    echo -e "id" > header.txt
    for file in !{name_file}
    do
         paste -d"\t" header.txt $file >temp.txt
         cat temp.txt > header.txt 
    done
    
    cat header.txt id_list.txt > all_tools_merge.matrix
    '''
}


/*
* Completion e-mail notification
*/
params.email = false

workflow.onComplete {

println print_cyan( workflow.success ? "Done!" : "Oops .. something went wrong" )
// Set up the e-mail variables
/*   def subject = "[nf-core/cirpipe] Successful: $workflow.runName"
if(!workflow.success){
    subject = "[nf-core/cirpipe] FAILED: $workflow.runName"
}
def email_fields = [:]
email_fields['version'] = workflow.manifest.version
email_fields['runName'] = custom_runName ?: workflow.runName
email_fields['success'] = workflow.success
email_fields['dateComplete'] = workflow.complete
email_fields['duration'] = workflow.duration
email_fields['exitStatus'] = workflow.exitStatus
email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
email_fields['errorReport'] = (workflow.errorReport ?: 'None')
email_fields['commandLine'] = workflow.commandLine
email_fields['projectDir'] = workflow.projectDir
email_fields['summary'] = summary
email_fields['summary']['Date Started'] = workflow.start
email_fields['summary']['Date Completed'] = workflow.complete
email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

// Render the TXT template
def engine = new groovy.text.GStringTemplateEngine()
def tf = new File("$baseDir/assets/email_template.txt")
def txt_template = engine.createTemplate(tf).make(email_fields)
def email_txt = txt_template.toString()

// Render the HTML template
def hf = new File("$baseDir/assets/email_template.html")
def html_template = engine.createTemplate(hf).make(email_fields)
def email_html = html_template.toString()

// Render the sendmail template
def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
def sf = new File("$baseDir/assets/sendmail_template.txt")
def sendmail_template = engine.createTemplate(sf).make(smail_fields)
def sendmail_html = sendmail_template.toString()

// Send the HTML e-mail
if (params.email) {
    try {
        if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
        // Try to send HTML e-mail using sendmail
        [ 'sendmail', '-t' ].execute() << sendmail_html
        log.info "[nf-core/cirpipe] Sent summary e-mail to $params.email (sendmail)"
    } catch (all) {
        // Catch failures and try with plaintext
        [ 'mail', '-s', subject, params.email ].execute() << email_txt
        log.info "[nf-core/cirpipe] Sent summary e-mail to $params.email (mail)"
    }
}

// Write summary e-mail HTML to a file
def output_d = new File( "${params.outdir}/Documentation/" )
if( !output_d.exists() ) {
    output_d.mkdirs()
}
def output_hf = new File( output_d, "pipeline_report.html" )
output_hf.withWriter { w -> w << email_html }
def output_tf = new File( output_d, "pipeline_report.txt" )
output_tf.withWriter { w -> w << email_txt }

log.info "[nf-core/cirpipe] Pipeline Complete"

}




def create_workflow_summary(summary) {

def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
yaml_file.text  = """
id: 'nf-core-cirpipe-summary'
description: " - this information is collected when the pipeline is started."
section_name: 'nf-core/cirpipe Workflow Summary'
section_href: 'https://github.com/nf-core/cirpipe'
plot_type: 'html'
data: |
    <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
    </dl>
""".stripIndent()

return yaml_file
*/

}


/*
* Parse software version numbers

process get_software_versions {

output:
file 'software_versions_mqc.yaml' into software_versions_yaml

script:
"""
echo $workflow.manifest.version > v_pipeline.txt
echo $workflow.nextflow.version > v_nextflow.txt
//multiqc --version > v_multiqc.txt
scrape_software_versions.py > software_versions_mqc.yaml
"""
}
*/





/*
* STEP 3 - Output Description HTML
process output_documentation {
tag "$prefix"
publishDir "${params.outdir}/Documentation", mode: 'copy'

input:
file output_docs

output:
file "results_description.html"

script:
"""
markdown_to_html.r $output_docs results_description.html
"""
}
*/



