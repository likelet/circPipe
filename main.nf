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

println(PATH)

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
     nf-core/circpipe v${workflow.manifest.version}
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
//summary['Run Name']     = custom_runName ?: workflow.runName
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

/*
def fastpoutdir = new File( "${params.outdir}/QC/" )
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

def autocircoutdir = new File( "${params.outdir}/pipeline_tools/pipeline_autocirc/" )
if( !autocircoutdir.exists() ) {
    autocircoutdir.mkdirs()
}

def tophatoutdir = new File( "${params.outdir}/pipeline_tools/pipeline_tophat/" )
if( !tophatoutdir.exists() ) {
    tophatoutdir.mkdirs()
}

def tophatcircexplorer2outdir = new File( "${params.outdir}/pipeline_tools/pipeline_circexplorer2_for_tophat/" )
if( !tophatcircexplorer2outdir.exists() ) {
    tophatcircexplorer2outdir.mkdirs()
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

def SEPoutdir = new File( "${params.outdir}/plot_separate/" )
if( !SEPoutdir.exists() ) {
    SEPoutdir.mkdirs()
}

def ALLoutdir = new File( "${params.outdir}/plot_merge/" )
if( !ALLoutdir.exists() ) {
    ALLoutdir.mkdirs()
}

def REPORToutdir = new File( "${params.outdir}/REPORT/" )
if( !REPORToutdir.exists() ) {
    REPORToutdir.mkdirs()
}
*/
otherTools = file(params.otherTools) //the other tools directory
if( !otherTools.exists() ) exit 1, print_red("Missing other tools directory: ${otherTools}")

if(params.mRNA){
    mRNA = file(params.mRNA) //the mRNA file
    if( !mRNA.exists() ) exit 1, print_red("Missing mRNA expression file: ${mRNA}")

}


/*
========================================================================================
                            check the index directory
========================================================================================
*/
if(params.starindex){
    starindex = Channel
            .fromPath(params.starindex)
            .ifEmpty { exit 1, "STAR index not found: ${params.starindex}" }
}

if(params.bowtie2index){
    bowtie2index = Channel
            .fromPath(params.bowtie2index)
            .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bowtie2index}" }

    bowtie2index_fc = Channel
            .fromPath(params.bowtie2index)
            .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bowtie2index}" }
}

if(params.bowtieindex){
    bowtieindex = Channel
            .fromPath(params.bowtieindex)
            .ifEmpty { exit 1, "Bowtie index not found: ${params.bowtieindex}" }
}

if(params.bwaindex){
    bwaindex = Channel
            .fromPath(params.bwaindex)
            .ifEmpty { exit 1, "BWA index not found: ${params.bwaindex}" }
}

if(params.segindex){
    segindex = Channel
            .fromPath(params.segindex)
            .ifEmpty { exit 1, "Segemehl index not found: ${params.segindex}" }
}

if(params.knifeindex){
    knifeindex = Channel
            .fromPath(params.knifeindex)
            .ifEmpty { exit 1, "KNIFE index not found: ${params.knifeindex}" }
}


/*
========================================================================================
                         the reference directory
========================================================================================
*/
refdir = file(params.refdir) //the reference genome directory
if( !refdir.exists() ) exit 1, print_red("Missing Reference Genome Directory: ${refdir}")

refmapsplice = file(params.refmapsplice) //the mapsplice reference genome directory
if( !refmapsplice.exists() ) exit 1, print_red("Missing Mapsplice Reference Genome Directory: ${refmapsplice}")

annotationfile = file(params.annotationfile) //the annotationfile
if( !annotationfile.exists() ) exit 1, print_red("Missing annotation file: ${annotationfile}")

genomefile = file(params.genomefile) //the genomefile
if( !genomefile.exists() ) exit 1, print_red("Missing genome file: ${genomefile}")

gtffile = file(params.gtffile) //the annotationfile-gtf-format
if( !gtffile.exists() ) exit 1, print_red("Missing gtf annotation file: ${gtffile}")

bedfile = file(params.bedfile) //the annotationfile-bed-format
if( !bedfile.exists() ) exit 1, print_red("Missing bed annotation file: ${bedfile}")


/*
========================================================================================
                         the environment directory
========================================================================================

condadir = file(params.condadir) //the python3 environment
if( !condadir.exists() ) exit 1, print_red("Missing python3 environment: ${condadir}")

conda2dir = file(params.conda2dir) //the python2 environment
if( !conda2dir.exists() ) exit 1, print_red("Missing python2 environment: ${conda2dir}")
*/

/*
========================================================================================
                         the tools directory
========================================================================================
*/
mapsdir = file(params.mapsdir) //the mapsplice directory
if( !mapsdir.exists() ) exit 1, print_red("Missing Mapsplice Directory: ${mapsdir}")

segdir = file(params.mapsdir) //the segemehl directory
if( !segdir.exists() ) exit 1, print_red("Missing Segemehl Directory: ${segdir}")

ciridir = file(params.ciridir)
if( !genomefile.exists() ) exit 1, print_red("Missing CIRI Directory: ${ciridir}")

find_circdir = file(params.find_circdir)
if( !find_circdir.exists() ) exit 1, print_red("Missing find_circ Directory: ${find_circdir}")

knifedir = file(params.knifedir)
if( !knifedir.exists() ) exit 1, print_red("Missing KNIFE Directory: ${knifedir}")



/*
========================================================================================
                         select the analysis tools
========================================================================================
*/
if( params.selectTools ==~ /.*1.*/ ){
    params.circexplorer2 = true
}else{
    params.circexplorer2 = false
}
if( params.selectTools ==~ /.*2.*/ ){
    params.ciri = true
}else{
    params.ciri = false
}
if( params.selectTools ==~ /.*3.*/ ){
    params.find_circ = true
}else{
    params.find_circ = false
}
if( params.selectTools ==~ /.*4.*/ ){
    params.mapsplice = true
}else{
    params.mapsplice = false
}
if( params.selectTools ==~ /.*5.*/ ){
    params.segemehl = true
}else{
    params.segemehl = false
}
if( params.selectTools ==~ /.*6.*/ ){
    params.knife = true
}else{
    params.knife = false
}


/*
========================================================================================
                         checking the design and compare file
========================================================================================
*/
//design file
if(params.designfile) {
    designfile = file(params.designfile)
    if( !designfile.exists() ) exit 1, print_red("Design file not found: ${params.designfile}")
}
//compare.txt
if(params.comparefile){
    comparefile = file(params.comparefile)
    if( !comparefile.exists() ) exit 1, print_red("Compare file not found: ${params.comparefile}")
}







/*
========================================================================================
                         showing the process and files
========================================================================================
*/
log.info print_cyan("""
========================================================================
    ________                          _______
   |  ____  |                        |  ___  |
   | |    |_|   _                    | |   | |   _
   | |         |_|                   | |   | |  |_|
   | |          _    _  __   _____   | |___| |   _    ______    _____
   | |         | |  | |/ /  |  ___|  |  _____|  | |  |  __  |  |  _  |
   | |         | |  |   /   | |      | |        | |  | |  | |  | |_| |
   | |     _   | |  |  /    | |      | |        | |  | |  | |  |  ___|
   | |____| |  | |  | |     | |___   | |        | |  | |__| |  | |___ 
   |________|  |_|  |_|     |_____|  |_|        |_|  |  ____|  |_____|
                                                     | |
                                                     | |
                                                     |_|

 =======================================================================
         """)
        .stripIndent()
log.info print_purple("============You are running cirPipe with the following parameters===============")
log.info print_purple("Checking parameters ...")
log.info "\n"
log.info print_yellow("========Reads types=======")
log.info print_yellow("singleEnd :                     ") + print_green(params.singleEnd)
log.info "\n"
log.info print_yellow("========Tools selected========")
log.info print_yellow("Circexplorer2 :                 ") + print_green(params.circexplorer2)
log.info print_yellow("Find_circ :                     ") + print_green(params.find_circ)
log.info print_yellow("Ciri :                          ") + print_green(params.ciri)
log.info print_yellow("Mapsplice :                     ") + print_green(params.mapsplice)
log.info print_yellow("Segemehl :                      ") + print_green(params.segemehl)
log.info "\n"
log.info print_yellow("=========Input files selected=========")
log.info print_yellow("Reads :                         ") + print_green(params.reads)
log.info print_yellow("Annotation file :               ") + print_green(params.annotationfile)
log.info print_yellow("Genome file :                   ") + print_green(params.genomefile)
log.info print_yellow("Gtf file :                      ") + print_green(params.gtffile)
log.info print_yellow("Bed file :                      ") + print_green(params.bedfile)
log.info "\n"
log.info print_yellow("==========Output files directory=========")
log.info print_yellow("Output directory :              ") + print_green(params.outdir)
log.info "\n"
log.info "\n"
log.info print_purple("Start running...")



/*
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
*/

/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */

Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { read_pairs_fastp }


/*
 * PREPROCESSING - Build STAR index
 */
if(!params.starindex){
    process makeSTARindex {
        publishDir "${params.outdir}/reference_genome", mode: 'copy', overwrite: true

        input:
        file genomefile
        file gtffile

        output:
        file "starindex" into starindex

        script:
        """
        mkdir starindex
        STAR \
            --runMode genomeGenerate \
            --runThreadN ${task.cpus} \
            --sjdbGTFfile ${gtffile} \
            --genomeDir starindex/ \
            --genomeFastaFiles ${genomefile} \
            --sjdbOverhang 149
        """
    }
}

/*
 * PREPROCESSING - Build BWA index
 */
if(!params.bwaindex){
    process makeBWAindex {
        publishDir "${params.outdir}/reference_genome", mode: 'copy', overwrite: true

        input:
        file genomefile


        output:
        file "bwaindex" into bwaindex

        script:
        """
        mkdir bwaindex
        cd ./bwaindex
        bwa \
            index ../${genomefile} \
            -p genome
        cd ../
        """
    }
}

/*
 * PREPROCESSING - Build Bowtie2 index
 */
if(!params.bowtie2index){
    process makeBowtie2index {
        publishDir "${params.outdir}/reference_genome", mode: 'copy', overwrite: true

        input:
        file genomefile


        output:
        file "bowtie2index" into bowtie2index
        file "bowtie2index" into bowtie2index_fc
        file "bowtie2index" into bowtie2_build_knife

        script:
        """
        mkdir bowtie2index
        cd ./bowtie2index
        bowtie2-build -f \
            ../${genomefile} \
            genome
        cd ../
        """
    }
}

/*
 * PREPROCESSING - Build Bowtie index
 */
if(!params.bowtieindex){
    process makeBowtieindex {
        publishDir "${params.outdir}/reference_genome", mode: 'copy', overwrite: true

        input:
        file genomefile


        output:
        file "bowtieindex" into bowtieindex
        file "bowtieindex" into bowtie_build_knife

        script:
        """
        mkdir bowtieindex
        cd ./bowtieindex
        bowtie-build \
            ../${genomefile} \
            genome
        cd ../
        """
    }
}

/*
 * PREPROCESSING - Build Segemehl index
 */
if(!params.segindex){
    process makeSegemehlindex {
        publishDir "${params.outdir}/reference_genome", mode: 'copy', overwrite: true

        input:
        file genomefile
        file segdir

        output:
        file "genome.idx" into segindex

        script:
        """
        ${segdir}/segemehl.x \
            -d ${genomefile} \
            -x genome.idx
        """
    }
}

/*
 * PREPROCESSING - Build KNIFE index
 */
if(!params.knifeindex){
    process makeKNIFEindex {

        input:
        file (bowtie_file) from bowtie_build_knife
        file (bowtie2_file) from bowtie2_build_knife
        file knifedir
        file gtffile

        output:
        file "KNIFE" into knife_use

        script:
        """
        mkdir KNIFE/circularRNApipeline_Standalone/denovo_scripts/index
        cp bowtieindex/* KNIFE/circularRNApipeline_Standalone/denovo_scripts/index
        mkdir KNIFE/circularRNApipeline_Standalone/index
        cp bowtie2index/* KNIFE/circularRNApipeline_Standalone/index
        cp ${gtffile} KNIFE/circularRNApipeline_Standalone/denovo_scripts
        """
    }
}




/*
========================================================================================
                       first step : run the fastp (QC tool)
========================================================================================
*/
process Fastp{
    tag "$pair_id"
    publishDir "${params.outdir}/QC", mode: 'copy', pattern: "*_fastpreport.html", overwrite: true

    input:
    set pair_id, file(query_file) from read_pairs_fastp

    output:
    set pair_id, file ('unzip_fastp_*') into fastpfiles_mapsplice,fastpfiles_bwa,fastpfiles_star,fastpfiles_segemehl,fastpfiles_tophat,fastpfiles_bowtie2
    file ('*.html') into fastp_for_waiting
    file ('*_fastp.json') into fastp_for_multiqc

    script:
    if(params.singleEnd){
        """
        fastp \
        -i ${query_file} \
        -o unzip_fastp_${pair_id}.fq \
        -h ${pair_id}_fastpreport.html \
        -j ${pair_id}_fastp.json
        """
    }else{
        """
        fastp \
        -i ${query_file[0]} \
        -I ${query_file[1]} \
        -o unzip_fastp_${pair_id}_1.fq \
        -O unzip_fastp_${pair_id}_2.fq \
        -h ${pair_id}_fastpreport.html \
        -j ${pair_id}_fastp.json 
        """
    }


}

fastp_for_waiting = fastp_for_waiting.first() //wait for finish this process first

/*
========================================================================================
                    run the multiqc (merge the results of fastp)
========================================================================================
*/
process Multiqc{
    publishDir "${params.outdir}/QC", mode: 'copy', pattern: "*.html", overwrite: true

    input:
    file (query_file) from fastp_for_multiqc.collect()

    output:
    file ('*.html') into multiqc_results

    script:
    """
    multiqc .
    """
}


/*
========================================================================================
                         the first tool : star - circexplorer2
                                      run the star
========================================================================================
*/
process Star{
    tag "$pair_id"
    publishDir "${params.outdir}/Alignment/STAR", mode: 'link', overwrite: true

    input:
    set pair_id, file(query_file) from fastpfiles_star
    file index from starindex.collect()

    output:
    set pair_id, file ('*.junction') into starfiles

    when:
    params.circexplorer2

    shell:
    if(params.singleEnd){
        """
        STAR \
        --runThreadN ${task.cpus} \
        --chimSegmentMin 10 \
        --genomeDir $index \
        --readFilesIn ${query_file} \
        --outFileNamePrefix star_${pair_id}_
        """
    }else{
        """
        STAR \
        --runThreadN ${task.cpus} \
        --chimSegmentMin 10 \
        --genomeDir $index \
        --readFilesIn ${query_file[0]} ${query_file[1]} \
        --outFileNamePrefix star_${pair_id}_
        """
    }

}

/*
========================================================================================
                         the first tool : star - circexplorer2
                                 run the circexplorer2
========================================================================================
*/
process Circexplorer2{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/CIRCexplorer2", mode: 'copy', overwrite: true

    input:
    set pair_id, file (query_file) from starfiles
    file annotationfile
    file genomefile

    output:
    set pair_id, file ('*known.txt') into circexplorer2files

    when:
    params.circexplorer2

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

/*
========================================================================================
                         the first tool : star - circexplorer2
                                 produce the bed6 file
========================================================================================
*/
process Circexplorer2_Bed{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/CIRCexplorer2", mode: 'copy', pattern: "*candidates.bed", overwrite: true

    input:
    set pair_id, file (query_file) from circexplorer2files
    file otherTools

    output:
    file ('*candidates.bed') into modify_circexplorer2
    val (pair_id) into modify_circexplorer2_id

    when:
    params.circexplorer2

    shell :
    '''
    grep circ !{query_file} \
    | grep -v chrM \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" "circexplorer2" "\t" $13 "\t" $6}' \
    > !{pair_id}_modify_circexplorer2.temp.bed
    
    python !{otherTools}/quchong.py !{pair_id}_modify_circexplorer2.temp.bed !{pair_id}_modify_circexplorer2.candidates.bed
    '''
}

/*
========================================================================================
                         the first tool : star - circexplorer2
                                  produce the matrix
========================================================================================
*/
process Circexplorer2_Matrix{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/CIRCexplorer2", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
    file (query_file) from modify_circexplorer2.collect()
    val (pair_id) from modify_circexplorer2_id.collect()
    file otherTools
    file designfile
    file gtffile

    output:
    file ('circexplorer2.txt') into merge_circexplorer2
    file ('*.matrix') into output_circexplorer2
    file ('name_circexplorer2.txt') into name_circexplorer2
    file ('*annote.txt') into de_circexplorer2
    file ('*.matrix') into plot_circexplorer2

    when:
    params.circexplorer2

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
        
    cat circexplorer2.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{otherTools}/bed1114.jar -i annotation.bed -o circexplorer2_ -gtf !{gtffile} -uniq

    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_circexplorer2.candidates.bed | wc -l`)) == 0 ];then
        python !{otherTools}/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{otherTools}/quchongsamples.py mergeconcatenate.bed ${line}_modify_circexplorer2.candidates.bed counts.txt
        fi
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
        echo -e "\\t${line}\\c" >> merge_header.txt
    done   
    
    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > circexplorer2_merge.matrix
    echo -e "circexplorer2" > name_circexplorer2.txt
    '''
}

/*
========================================================================================
                          the first tool : star - circexplorer2
                                 Differential Expression
========================================================================================
*/
process Circexplorer2_DE{
    publishDir "${params.outdir}/DE_Analysis/Circexplorer2", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (anno_file) from de_circexplorer2
    file otherTools
    file designfile
    file comparefile
    file (matrix_file) from plot_circexplorer2

    output:
    file ('*') into end_circexplorer2

    when:
    params.circexplorer2

    shell:
    '''
    Rscript !{otherTools}/edgeR_circ.R !{otherTools}/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}


/*
========================================================================================
                              the second tool : bwa - ciri
                                      run the bwa
========================================================================================
*/
process Bwa{
    tag "$pair_id"
    publishDir "${params.outdir}//Alignment/BWA", mode: 'link', overwrite: true

    input:
    set pair_id, file (query_file) from fastpfiles_bwa
    file index from bwaindex.collect()

    output:
    set pair_id, file ('*.sam') into bwafiles

    when:
    params.ciri

    shell:
    if(params.singleEnd){
        """
        bwa \
        mem -t ${task.cpus} \
        -k 15 \
        -T 19  -M -R \
        "@RG\\tID:fastp_${pair_id}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${pair_id}" \
        ${index}/genome \
        ${query_file} \
        > bwa_${pair_id}.mem.sam
        """
    }else{
        """
        bwa \
        mem -t ${task.cpus} \
        -T 19 -M -R \
        "@RG\\tID:fastp_${pair_id}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${pair_id}" \
        ${index}/genome \
        ${query_file[0]} ${query_file[1]} \
        > bwa_${pair_id}.mem.sam
        """
    }

}

/*
========================================================================================
                              the second tool : bwa - ciri
                                      run the ciri
========================================================================================
*/
process Ciri{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/CIRI", mode: 'copy', overwrite: true

    input:
    set pair_id, file (query_file) from bwafiles
    file gtffile
    file genomefile
    file ciridir

    output:
    set pair_id, file ('*.txt') into cirifiles

    when:
    params.ciri

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

/*
========================================================================================
                              the second tool : bwa - ciri
                                  produce the bed6 file
========================================================================================
*/
process Ciri_Bed{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/CIRI", mode: 'copy', pattern: "*candidates.bed", overwrite: true

    input:
    set pair_id, file (query_file) from cirifiles
    file otherTools

    output:
    file ('*candidates.bed') into modify_ciri_file
    val (pair_id) into modify_ciri_id

    when:
    params.ciri

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

/*
========================================================================================
                              the second tool : bwa - ciri
                                  produce the matrix
========================================================================================
*/
process Ciri_Matrix{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/CIRI", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
    file (query_file) from modify_ciri_file.collect()
    val (pair_id) from modify_ciri_id.collect()
    file otherTools
    file designfile
    file gtffile

    output:
    file ('ciri.txt') into merge_ciri
    file ('*.matrix') into output_ciri
    file ('name_ciri.txt') into name_ciri
    file ('*annote.txt') into de_ciri
    file ('*.matrix') into plot_ciri

    when:
    params.ciri

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
        
    cat ciri.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{otherTools}/bed1114.jar -i annotation.bed -o ciri_ -gtf !{gtffile} -uniq

    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_ciri.candidates.bed | wc -l`)) == 0 ];then
        python !{otherTools}/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{otherTools}/quchongsamples.py mergeconcatenate.bed ${line}_modify_ciri.candidates.bed counts.txt
        fi
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
        echo -e "\\t${line}\\c" >> merge_header.txt
    done   
       
    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > ciri_merge.matrix
    echo -e "ciri" > name_ciri.txt
    '''
}

/*
========================================================================================
                               the second tool : bwa - ciri
                                 Differential Expression
========================================================================================
*/
process Ciri_DE{
    publishDir "${params.outdir}/DE_Analysis/CIRI", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (anno_file) from de_ciri
    file otherTools
    file designfile
    file comparefile
    file (matrix_file) from plot_ciri

    output:
    file ('*') into end_ciri

    when:
    params.ciri

    shell:
    '''
    Rscript !{otherTools}/edgeR_circ.R !{otherTools}/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}


/*
========================================================================================
                              the third tool : mapsplice
                                  run the mapsplice
========================================================================================
*/
process Mapsplice{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Mapsplice", mode: 'copy', overwrite: true

    input:
    set pair_id, file (query_file) from fastpfiles_mapsplice
    file mapsdir
    file gtffile
    file refmapsplice
    file outdir
    file index from bowtieindex.collect()

    output:
    set pair_id, file('*') into mapsplicefiles


    when:
    params.mapsplice

    shell:
    if(params.singleEnd){
        """
        source activate tools_in_python3
        python ${mapsdir}/mapsplice.py \
        -p ${task.cpus} \
        -k 1 \
        --fusion-non-canonical \
        --qual-scale phred33 \
        --non-canonical-double-anchor \
        --min-fusion-distance 200 \
        -x ${index}/genome \
        --gene-gtf ${gtffile} \
        -c ${refmapsplice} \
        -1 ${query_file} \
        -o output_mapsplice_${pair_id} 2 \
        > ${pair_id}_mapsplice.log
        source deactivate      
        """
    }else{
        """
        source activate tools_in_python3
        python ${mapsdir}/mapsplice.py \
        -p ${task.cpus} \
        -k 1 \
        --fusion-non-canonical \
        --non-canonical-double-anchor \
        --min-fusion-distance 200 \
        -x ${index}/genome \
        --gene-gtf ${gtffile} \
        -c ${refmapsplice} \
        -1 ${query_file[0]} \
        -2 ${query_file[1]} \
        -o output_mapsplice_${pair_id} 2 \
        > ${pair_id}_mapsplice.log
        source deactivate
        """
    }

}

/*
========================================================================================
                              the third tool : mapsplice
                                 produce the bed6 file
========================================================================================
*/
process Mapsplice_Bed{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Mapsplice", mode: 'copy', pattern: "*candidates.bed", overwrite: true

    input:
    set pair_id, file (query_file) from mapsplicefiles
    file outdir
    file otherTools

    output:
    file ('*candidates.bed') into modify_mapsplice
    val (pair_id) into modify_mapsplice_id

    when:
    params.mapsplice

    shell :
    '''
    if [ $((`cat output_mapsplice_!{pair_id}/circular_RNAs.txt | wc -l`)) == 0 ];then
    touch !{pair_id}_modify_mapsplice.candidates.bed
    else
    cat output_mapsplice_!{pair_id}/circular_RNAs.txt \
    | awk '{print $6}' \
    | sed -e 's/.//' \
    > !{pair_id}_mapsplice_temp1.bed

    cat output_mapsplice_!{pair_id}/circular_RNAs.txt \
    | awk '{print $1}' \
    | awk -F"~" '{print $2}' \
    > !{pair_id}_mapsplice_temp.bed

    paste !{pair_id}_mapsplice_temp.bed !{pair_id}_mapsplice_temp1.bed output_mapsplice_!{pair_id}/circular_RNAs.txt \
    | grep -v chrM \
    | awk '{if($2=="-") print $1 "\t" $4 "\t" $5 "\t" "mapsplice" "\t" $7 "\t" $2 ; else print $1 "\t" $5 "\t" $4 "\t" "mapsplice" "\t" $7 "\t" $2 }' \
    > !{pair_id}_modify_mapsplice.temp.bed
    
    python !{otherTools}/quchong.py !{pair_id}_modify_mapsplice.temp.bed !{pair_id}_modify_mapsplice.candidates.bed
    fi
    '''
}

/*
========================================================================================
                              the third tool : mapsplice
                                  produce the matrix
========================================================================================
*/
process Mapsplice_Matrix{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Mapsplice", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
    file (query_file) from modify_mapsplice.collect()
    val (pair_id) from modify_mapsplice_id.collect()
    file otherTools
    file designfile
    file gtffile

    output:
    file ('mapsplice.txt') into merge_mapsplice
    file ('*.matrix') into output_mapsplice
    file ('name_mapsplice.txt') into name_mapsplice
    file ('*annote.txt') into de_mapsplice
    file ('*.matrix') into plot_mapsplice

    when:
    params.mapsplice

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
        
    cat mapsplice.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{otherTools}/bed1114.jar -i annotation.bed -o mapsplice_ -gtf !{gtffile} -uniq
    
    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_mapsplice.candidates.bed | wc -l`)) == 0 ];then
        python !{otherTools}/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{otherTools}/quchongsamples.py mergeconcatenate.bed ${line}_modify_mapsplice.candidates.bed counts.txt
        fi
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
        echo -e "\\t${line}\\c" >> merge_header.txt
    done   

    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > mapsplice_merge.matrix
    echo -e "mapsplice" > name_mapsplice.txt
    '''
}

/*
========================================================================================
                                the third tool : mapsplice
                                 Differential Expression
========================================================================================
*/
process Mapsplice_DE{
    publishDir "${params.outdir}/DE_Analysis/Mapsplice", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (anno_file) from de_mapsplice
    file otherTools
    file designfile
    file comparefile
    file (matrix_file) from plot_mapsplice

    output:
    file ('*') into end_mapsplice

    when:
    params.mapsplice

    shell:
    '''
    Rscript !{otherTools}/edgeR_circ.R !{otherTools}/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}


/*
========================================================================================
                              the fourth tool : segemehl
                                   run the segemehl
========================================================================================
*/
process Segemehl{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Segemehl", mode: 'copy', overwrite: true

    input:
    set pair_id, file (query_file) from fastpfiles_segemehl
    file segdir
    file genomefile
    file index from segindex.collect()

    output:
    set pair_id, file ('*splicesites.bed') into segemehlfiles


    when:
    params.segemehl

    shell:
    if(params.singleEnd){
        """
        source activate tools_in_python3
        ${segdir}/segemehl.x \
        -d ${genomefile} \
        -i ${index} \
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
        source deactivate
        """
    }else{
        """
        source activate tools_in_python3
        ${segdir}/segemehl.x \
        -d ${genomefile} \
        -i ${index} \
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
        source deactivate
        """
    }

}

/*
========================================================================================
                              the fourth tool : segemehl
                                 produce the bed6 file
========================================================================================
*/
process Segemehl_Bed{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Segemehl", mode: 'copy', pattern:"*candidates.bed", overwrite: true

    input:
    set pair_id , file ( query_file ) from segemehlfiles
    file otherTools

    output:
    file ('*candidates.bed') into modify_segemehl
    val (pair_id) into modify_segemehl_id

    when:
    params.segemehl

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

/*
========================================================================================
                              the fourth tool : segemehl
                                  produce the matrix
========================================================================================
*/
process Segemehl_Matrix{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Segemehl", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
    file (query_file) from modify_segemehl.collect()
    val (pair_id) from modify_segemehl_id.collect()
    file otherTools
    file designfile
    file gtffile

    output:
    file ('segemehl.txt') into merge_segemehl
    file ('*.matrix') into output_segemehl
    file ('*annote.txt') into de_segemehl
    file ('*.matrix') into plot_segemehl
    file ('name_segemehl.txt') into name_segemehl

    when:
    params.segemehl

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
    
    cat segemehl.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{otherTools}/bed1114.jar -i annotation.bed -o segemehl_ -gtf !{gtffile} -uniq

    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_segemehl.candidates.bed | wc -l`)) == 0 ];then
        python !{otherTools}/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{otherTools}/quchongsamples.py mergeconcatenate.bed ${line}_modify_segemehl.candidates.bed counts.txt
        fi
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
        echo -e "\\t${line}\\c" >> merge_header.txt
    done   
       
    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > segemehl_merge.matrix
    echo -e "segemehl" > name_segemehl.txt
    '''
}

/*
========================================================================================
                               the fourth tool : segemehl
                                Differential Expression
========================================================================================
*/
process Segemehl_DE{
    publishDir "${params.outdir}/DE_Analysis/Segemehl", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (anno_file) from de_segemehl
    file otherTools
    file designfile
    file comparefile
    file (matrix_file) from plot_segemehl

    output:
    file ('*') into end_segemehl

    when:
    params.segemehl

    shell:
    '''
    Rscript !{otherTools}/edgeR_circ.R !{otherTools}/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}


/*
========================================================================================
                          the fifth tool : bowtie2 - find_circ
                                   run the bowtie2
========================================================================================
*/
process Bowtie2{
    tag "$pair_id"
    publishDir "${params.outdir}/Alignment/Bowtie2", mode: 'link', overwrite: true

    input:
    set pair_id, file (query_file) from fastpfiles_bowtie2
    file index from bowtie2index.collect()

    output:
    set pair_id, file ('bowtie2_unmapped_*') into bowtie2files
    set pair_id, file ('bowtie2_unmapped_*') into bowtie2files_for_autocirc

    when:
    params.find_circ

    shell:
    if(params.singleEnd){
        """
        bowtie2 \
        -p ${task.cpus} \
        --very-sensitive \
        --score-min=C,-15,0 \
        --mm \
        -x ${index}/genome \
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
        -x ${index}/genome \
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

/*
========================================================================================
                          the fifth tool : bowtie2 - find_circ
                                   run the find_circ
========================================================================================
*/
process Find_circ{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Find_circ", mode: 'copy', overwrite: true

    input:
    set pair_id, file (query_file) from bowtie2files
    file genomefile
    file find_circdir
    file index from bowtie2index_fc.collect()

    output:
    set pair_id, file ('*splice_sites.bed') into find_circfiles


    when:
    params.find_circ

    shell:
    """     
    source activate tools_in_python2
    python ${find_circdir}/unmapped2anchors.py ${query_file} \
    | gzip \
    > find_circ_${pair_id}_anchors.qfa.gz

    bowtie2 \
    -p ${task.cpus} \
    --reorder \
    --mm \
    --score-min=C,-15,0 \
    -q \
    -x ${index}/genome \
    -U find_circ_${pair_id}_anchors.qfa.gz \
    | python ${find_circdir}/find_circ.py \
    -G ${genomefile} \
    -p ${pair_id}_ \
    -s find_circ_${pair_id}_stats.sites.log \
    -n find_circ \
    -R find_circ_${pair_id}_spliced_reads.fa \
    > find_circ_${pair_id}_splice_sites.bed   

    source deactivate
    """
}

/*
========================================================================================
                          the fifth tool : bowtie2 - find_circ
                                 produce the bed6 file
========================================================================================
*/
process Find_circ_Bed{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Find_circ", mode: 'copy', pattern: "*candidates.bed", overwrite: true

    input:
    set pair_id, file (query_file) from find_circfiles
    file otherTools

    output:
    file ('*candidates.bed') into modify_find_circfiles
    val (pair_id) into modify_find_circ_id

    when:
    params.find_circ

    shell :
    '''
    grep CIRCULAR !{query_file} \
    | grep -v chrM \
    | grep UNAMBIGUOUS_BP \
    | grep ANCHOR_UNIQUE \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $5 "\t" $6}' \
    | sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 \
    > !{pair_id}_modify_find_circ.temp.bed
    
    python !{otherTools}/quchong.py !{pair_id}_modify_find_circ.temp.bed !{pair_id}_modify_find_circ.candidates.bed
    '''
}

/*
========================================================================================
                          the fifth tool : bowtie2 - find_circ
                                 produce the matrix
========================================================================================
*/

process Find_circ_Matrix{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/Find_circ", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
    file (query_file) from modify_find_circfiles.collect()
    val (pair_id) from modify_find_circ_id.collect()
    file otherTools
    file designfile
    file gtffile

    output:
    file ('find_circ.txt') into merge_find_circ
    file ('*annote.txt') into de_find_circ
    file ('*.matrix') into output_find_circ
    file ('*.matrix') into plot_find_circ
    file ('name_find_circ.txt') into name_find_circ
    file ('*annote.txt') into cor_find_circ
    file ('*.matrix') into plot_find_circ_cor


    when:
    params.find_circ

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
    
    cat find_circ.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{otherTools}/bed1114.jar -i annotation.bed -o find_circ_ -gtf !{gtffile} -uniq

    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_find_circ.candidates.bed | wc -l`)) == 0 ];then
        python !{otherTools}/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{otherTools}/quchongsamples.py mergeconcatenate.bed ${line}_modify_find_circ.candidates.bed counts.txt
        fi
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
        echo -e "\\t${line}\\c" >> merge_header.txt
    done   
    
    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > find_circ_merge.matrix
    echo -e "find_circ" > name_find_circ.txt
    '''
}

/*
========================================================================================
                          the fifth tool : bowtie2 - find_circ
                                 Differential Expression
========================================================================================
*/
process Find_circ_DE{
    publishDir "${params.outdir}/DE_Analysis/Find_circ", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (anno_file) from de_find_circ
    file otherTools
    file designfile
    file comparefile
    file (matrix_file) from plot_find_circ

    output:
    file ('*') into end_find_circ

    when:
    params.find_circ

    shell:
    '''
    Rscript !{otherTools}/edgeR_circ.R !{otherTools}/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}

/*
========================================================================================
                          the fifth tool : bowtie2 - find_circ
                                        Correlation
========================================================================================
*/
process Find_circ_Cor{
    publishDir "${params.outdir}/Corrrelation_Analysis/Find_circ", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (matrix_file) from plot_find_circ_cor
    file (anno_file) from cor_find_circ
    file mRNA
    file otherTools

    when:
    params.mRNA && params.find_circ

    output:
    file ('*') into cor_plot

    shell:
    '''
    Rscript !{otherTools}/correlation.R !{otherTools}/R_function.R !{mRNA} !{matrix_file} !{anno_file}
    '''
}


/*
========================================================================================
                                after running the tools
                calculate the results by different tools, combine matrix
========================================================================================
*/
process Tools_Merge{
    publishDir "${params.outdir}/Combination_Matrix", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
    file (query_file) from merge_find_circ.concat( merge_circexplorer2, merge_ciri, merge_mapsplice, merge_segemehl ).collect()
    file (name_file) from name_find_circ.concat( name_circexplorer2, name_ciri, name_mapsplice, name_segemehl ).collect()
    file (matrix_file) from output_find_circ.concat( output_circexplorer2, output_ciri, output_mapsplice, output_segemehl ).collect()
    file otherTools

    output:
    file ('all_tools_merge.matrix') into tools_merge
    file ('for_annotation.bed') into bed_for_annotation
    file ('final.matrix') into matrix_for_circos

    when:
    params.merge

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
    
    for file in !{matrix_file}
    do
        cat $file | awk 'NR==1' > sample_id.txt
        cat $file > temp.txt
        sed -i '1d' temp.txt
        cat temp.txt >> total_matrix.txt
    done
    
    Rscript !{otherTools}/changematrix.R total_matrix.txt change_reads.txt
    
    python !{otherTools}/finalmerge.py change_reads.txt newmatrix.txt for_annotation.bed

    cat sample_id.txt newmatrix.txt > final.matrix
    '''
}

/*
========================================================================================
                                after running the tools
                                        annoation
========================================================================================
*/
process Annotation_Merge{
    publishDir "${params.outdir}/Annotation", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (bed_file) from bed_for_annotation
    file (query_file) from matrix_for_circos
    file otherTools
    file gtffile

    when:
    params.merge

    output:
    file ('*') into annotation_plot

    shell:
    '''
    java -jar !{otherTools}/bed1114.jar -i !{bed_file} -o merge_ -gtf !{gtffile} -uniq 
    Rscript !{otherTools}/circos.R !{query_file}
    Rscript !{otherTools}/circRNA_feature.R !{otherTools}/R_function.R merge_for_annotation_annote.txt
    '''
}

process Venn{
    publishDir "${params.outdir}/Annotation", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (matrix_file) from tools_merge


    when:
    params.merge

    output:
    file ('*') into annotation_plot

    shell:
    '''
    Rscript !{otherTools}/venn.R !{matrix_file} venn.png
    '''
}



/*
* Completion e-mail notification
*/

emailaddress = params.email

/*
params.plaintext_email = false
if(params.email) {
    summary['E-mail Address'] = params.email
}
*/

workflow.onComplete {

println print_cyan( workflow.success ? "Done!" : "Oops .. something went wrong" )

    def msg = """\
<html>
<head>
  <meta charset="utf-8">
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="circpipe: cirRNA analysis pipe">
  <title>circpipe Pipeline Report</title>
</head>
<body>
<div style="text-align:center; font-family: Helvetica, Arial, sans-serif; padding: 30px; max-width: 800px; margin: 0 auto;">
    <h1> Pipeline execution summary </h1>
    <h2> ---------------------------------------------- </h2>
    <div style = "text-align:center; color: #3c763d; background-color: #dff0d8; border-color: #d6e9c6; padding: 15px; margin-bottom: 20px; border: 1px solid transparent; border-radius: 4px;" >
        <h3 style = "margin-top:0; color: inherit;" ><strong>CircPipe</strong></h3>
        <p>Completed at : <strong>${workflow.complete}</strong></p>
        <p>Duration : <strong>${workflow.duration}</strong></p>
        <p>Success : <strong>${workflow.success}</strong></p>
        <p>Exit status : <strong>${workflow.exitStatus}</strong></p>
    </div>
    <div style="text-align:center; color: #a94442; background-color: #f2dede; border-color: #ebccd1; padding: 15px; margin-bottom: 20px; border: 1px solid transparent; border-radius: 4px;">
        <p>Error report : </p>
        <pre style="white-space: pre-wrap; overflow: visible; margin-bottom: 0;">${workflow.errorReport}</pre>
    </div>
    <p>The command used to launch the workflow was as follows : </p>      
    <pre style="white-space: pre-wrap; overflow: visible; background-color: #ededed; padding: 15px; border-radius: 4px; margin-bottom:0px;">${workflow.commandLine}</pre>
    <h3> Tools selected : </h3>
    <table style="width:100%; max-width:100%; border-spacing: 0; border-collapse: collapse; border:0; margin-bottom: 30px;">
    <tbody style="border-bottom: 1px solid #ddd;">
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Circexplorer2 : ${params.circexplorer2} </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:cneter; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Find_circ : ${params.find_circ} </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Ciri : ${params.ciri} </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Mapsplice : ${params.mapsplice} </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Segemehl : ${params.segemehl} </pre></td>
    </tr>
    </tbody>
    </table>

    <p> likelet/circPipe </p>
    <p><a href="https://github.com/likelet/cirPipe">https://github.com/likelet/circPipe</a></p>
</div>
</body>
</html>
        """
            .stripIndent()

    sendMail(to: emailaddress,
            subject: 'Breaking News in CircPipe Mission!',
            body: msg)
}





// Set up the e-mail variables
/*   def subject = "[nf-core/cirpipe] Successful: $workflow.runName"
if(!workflow.success){
    subject = "[nf-core/cirpipe] FAILED: $workflow.runName"
}
def email_fields = [:]
email_fields['version'] = workflow.manifest.version
//email_fields['runName'] = custom_runName ?: workflow.runName
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
*/



/*

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


}
*/

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

workflow.onError {

    println print_cyan(workflow.success ? "Done!" : "Oops .. something went wrong")

    def msg = """\
<html>
<head>
  <meta charset="utf-8">
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="circpipe: cirRNA analysis pipe">
  <title>circpipe Pipeline Report</title>
</head>
<body>
<div style="text-align:center; font-family: Helvetica, Arial, sans-serif; padding: 30px; max-width: 800px; margin: 0 auto;">
    <h1> Pipeline execution summary </h1>
    <h2> ---------------------------------------------- </h2>
    <div style = "text-align:center; color: #3c763d; background-color: #dff0d8; border-color: #d6e9c6; padding: 15px; margin-bottom: 20px; border: 1px solid transparent; border-radius: 4px;" >
        <h3 style = "margin-top:0; color: inherit;" ><strong>CircPipe</strong></h3>
        <p>Completed at : <strong>${workflow.complete}</strong></p>
        <p>Duration : <strong>${workflow.duration}</strong></p>
        <p>Success : <strong>${workflow.success}</strong></p>
        <p>Exit status : <strong>${workflow.exitStatus}</strong></p>
    </div>
    <div style="text-align:center; color: #a94442; background-color: #f2dede; border-color: #ebccd1; padding: 15px; margin-bottom: 20px; border: 1px solid transparent; border-radius: 4px;">
        <p>Error report : </p>
        <pre style="white-space: pre-wrap; overflow: visible; margin-bottom: 0;">${workflow.errorMessage}</pre>
    </div>
    <p>The command used to launch the workflow was as follows : </p>      
    <pre style="white-space: pre-wrap; overflow: visible; background-color: #ededed; padding: 15px; border-radius: 4px; margin-bottom:0px;">${
        workflow.commandLine
    }</pre>
    <h3> Tools selected : </h3>
    <table style="width:100%; max-width:100%; border-spacing: 0; border-collapse: collapse; border:0; margin-bottom: 30px;">
    <tbody style="border-bottom: 1px solid #ddd;">
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Circexplorer2 : ${
        params.circexplorer2
    } </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:cneter; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Find_circ : ${
        params.find_circ
    } </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Ciri : ${
        params.ciri
    } </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Mapsplice : ${
        params.mapsplice
    } </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Segemehl : ${
        params.segemehl
    } </pre></td>
    </tr>
    </tbody>
    </table>

    <p> likelet/circPipe </p>
    <p><a href="https://github.com/likelet/cirPipe">https://github.com/likelet/circPipe</a></p>
</div>
</body>
</html>
        """
            .stripIndent()

    sendMail(to: '513848731@qq.com',
            subject: 'Breaking News in CircPipe Mission!',
            body: msg)
}

