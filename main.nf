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

// Configurable variables
params.name = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
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

/*
 * Create a channel for input read files
 */
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
summary['Fasta Ref']    = params.fasta
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




/*
 * Checking the input files
 * Adding input files error exceptions Here
 */
outdir = file(params.outdir) //the output directory
if( !outdir.exists() ) exit 1, print_red("Missing output directory: ${outdir}")

if ( params.star ){

    starindex = file(params.starindex) //the index directory
    if( !starindex.exists() ) exit 1, print_red("Missing star index directory: ${starindex}")

}

if ( params.mapsplice ){

    condadir = file(params.condadir) //the python3 environment
    if( !condadir.exists() ) exit 1, print_red("Missing python3 environment: ${condadir}")

    gtffile = file(params.gtffile) //the annotationfile
    if( !gtffile.exists() ) exit 1, print_red("Missing annotation file: ${gtffile}")

    mapsdir = file(params.mapsdir) //the mapsplice directory
    if( !mapsdir.exists() ) exit 1, print_red("Missing Mapsplice Directory: ${mapsdir}")

    refdir = file(params.refdir) //the reference genome directory
    if( !refdir.exists() ) exit 1, print_red("Missing Reference Genome Directory: ${refdir}")

}

if ( params.segemehl ){

    condadir = file(params.condadir) //the python3 environment
    if( !condadir.exists() ) exit 1, print_red("Missing python3 environment: ${condadir}")

    genomefile = file(params.genomefile) //the genomefile
    if( !genomefile.exists() ) exit 1, print_red("Missing genome file: ${genomefile}")

    segdir = file(params.mapsdir) //the segemehl directory
    if( !segdir.exists() ) exit 1, print_red("Missing Segemehl Directory: ${segdir}")

    segindex = file(params.segindex) //the segemehl index file
    if( !segindex.exists() ) exit 1, print_red("Missing Segemehl index file: ${segindex}")

}

if ( params.circexplorer2 ){

    if ( params.star ){

        annotationfile = file(params.annotationfile) //the annotationfile
        if( !annotationfile.exists() ) exit 1, print_red("Missing annotation file: ${annotationfile}")

        genomefile = file(params.genomefile) //the genomefile
        if( !genomefile.exists() ) exit 1, print_red("Missing genome file: ${genomefile}")

    }else{

        exit 1, print_yellow( "Please run star before circexplorer2" )

    }


}

if ( params.ciri ){

    if ( params.bwa ){

        gtffile = file(params.gtffile) //the annotationfile
        if( !gtffile.exists() ) exit 1, print_red("Missing annotation file: ${gtffile}")

        genomefile = file(params.genomefile) //the genomefile
        if( !genomefile.exists() ) exit 1, print_red("Missing genome file: ${genomefile}")

        ciridir = file(params.ciridir)
        if( !genomefile.exists() ) exit 1, print_red("Missing CIRI Directory: ${ciridir}")

    }else{

        exit 1, print_yellow( "Please run bwa before ciri" )

    }

}

if ( params.find_circ ){

    if ( params.bowtie2 ){

        conda2dir = file(params.conda2dir) //the python2 environment
        if( !conda2dir.exists() ) exit 1, print_red("Missing python2 environment: ${conda2dir}")

        genomefile = file(params.genomefile) //the genomefile
        if( !genomefile.exists() ) exit 1, print_red("Missing genome file: ${genomefile}")

        find_circdir = file(params.find_circdir)
        if( !find_circdir.exists() ) exit 1, print_red("Missing find_circ Directory: ${find_circdir}")

    }else{

        exit 1, print_yellow( "Please run bowtie2 before find_circ" )

    }

}


//showing the process and files
log.info print_purple("""\
         c i r P i p e   P I P E L I N E
         =============================


         Reads types :
         singleEnd : ${params.singleEnd}
         
         Tools selected :
         star : ${params.star}
         bwa : ${params.bwa}
         bowtie2 = ${params.bowtie2}
         mapsplice = ${params.mapsplice}
         segemehl = ${params.segemehl}
         find_circ = ${params.find_circ}
         circexplorer2 = ${params.circexplorer2}
         ciri = ${params.ciri}

         Input files selected :
         reads : ${params.reads}
         annotation file : ${params.annotationfile}
         genome file : ${params.genomefile}
         gtf file : ${params.gtffile}

         Index selected :
         segemehl index : ${params.segindex}

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
Channel
        .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { read_pairs_fastp }


//run the fastp
process run_fastp{
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy', pattern: "*.html", overwrite: true

    maxForks fork_number

    input:
    set pair_id, file(query_file) from read_pairs_fastp

    output:
    set pair_id, file ('fastp_*') into fastpfiles_star
    set pair_id, file ('fastp_*') into fastpfiles_bwa
    set pair_id, file ('unzip_fastp_*') into fastpfiles_mapsplice
    set pair_id, file ('unzip_fastp_*') into fastpfiles_segemehl
    set pair_id, file ('fastp_*') into fastpfiles_bowtie2
    file ('*.html') into fastqc_for_waiting


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
    -h ${pair_id}.html
    """

}

fastqc_for_waiting = fastqc_for_waiting.first() //wait for finish this process first


//start mapping
if ( params.star ){

    process run_star{
        tag "$pair_id"
        publishDir params.outdir, mode: 'link', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file(query_file) from fastpfiles_star
        file starindex

        output:
        set pair_id, file ('*.junction') into starfiles

        shell:
        star_threads = idv_cpu - 1
        """
        STAR \
    	--runThreadN ${star_threads} \
    	--chimSegmentMin 10 \
    	--genomeDir ${starindex} \
    	--readFilesCommand zcat \
    	--readFilesIn ${query_file[0]} ${query_file[1]} \
    	--outFileNamePrefix star_${pair_id}_
        """
    }

}

if ( params.bwa ){

    process run_bwa{
        tag "$pair_id"
        publishDir params.outdir, mode: 'link', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from fastpfiles_bwa

        output:
        set pair_id, file ('*.sam') into bwafiles

        shell:
        bwa_threads = idv_cpu - 1
        """
        bwa \
    	mem -t ${bwa_threads} \
    	-T 19 -M -R \
    	"@RG\\tID:fastp_${pair_id}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${pair_id}" \
    	/home/wqj/test/bwaindex/genome \
    	${query_file[0]} ${query_file[1]} \
	    > bwa_${pair_id}.mem.sam
        """
    }

}

if ( params.mapsplice ){

    process run_mapsplice{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from fastpfiles_mapsplice
        file mapsdir
        file gtffile
        file refdir
        file outdir

        output:
        set pair_id into mapsplicefiles

        conda params.condadir

        shell:
        mapsplice_threads = idv_cpu - 1
        """
	    python ${mapsdir}/mapsplice.py \
	    -p ${mapsplice_threads} \
	    -k 1 \
	    --fusion-non-canonical \
	    --non-canonical-double-anchor \
	    --min-fusion-distance 200 \
	    -x /home/wqj/test/bowtieindex/chrX \
	    --gene-gtf ${gtffile} \
	    -c ${refdir} \
	    -1 ${query_file[0]} \
	    -2 ${query_file[1]} \
	    -o ${outdir}/output_mapsplice_${pair_id}
	    """
    }

    process run_modify_mapsplice{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', pattern: "*candidates.bed", overwrite: true

        maxForks fork_number

        input:
        set pair_id from mapsplicefiles
        file outdir

        output:
        set pair_id, file ('*candidates.bed') into modify_mapsplice

        shell :
        '''
        cat !{outdir}/output_mapsplice_!{pair_id}/circular_RNAs.txt \
	    | awk '{print $6}' \
	    | sed -e 's/.//' \
	    > !{pair_id}_mapsplice_temp1.bed

	    cat !{outdir}/output_mapsplice_!{pair_id}/circular_RNAs.txt \
	    | awk '{print $1}' \
	    | awk -F"~" '{print $2}' \
	    > !{pair_id}_mapsplice_temp.bed

	    paste !{pair_id}_mapsplice_temp.bed !{pair_id}_mapsplice_temp1.bed !{outdir}/output_mapsplice_!{pair_id}/circular_RNAs.txt \
	    | grep -v chrM \
	    | awk '{if($2=="-") print $1 "\t" $4 "\t" $5 "\t" "mapsplice" "\t" $7 "\t" $2 ; else print $1 "\t" $5 "\t" $4 "\t" "mapsplice" "\t" $7 "\t" $2 }' \
	    > !{pair_id}_modify_mapsplice.candidates.bed
        '''
    }

}

if ( params.segemehl ){

    process run_segemehl{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from fastpfiles_segemehl
        file segdir
        file genomefile
        file segindex

        output:
        set pair_id, file ('*splicesites.bed') into segemehlfiles

        conda params.condadir

        shell:
        segemehl_threads = idv_cpu - 1
        """
	    ${segdir}/segemehl.x \
	    -d ${genomefile} \
	    -i ${segindex} \
	    -q ${query_file[0]} \
	    -p ${query_file[1]} \
	    -t ${segemehl_threads} \
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

    process run_modify_segemehl{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', pattern:"*candidates.bed", overwrite: true

        maxForks fork_number

        input:
        set pair_id , file ( query_file ) from segemehlfiles

        output:
        set pair_id, file ('*candidates.bed') into modify_segemehl

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
	    > !{pair_id}_modify_segemehl.candidates.bed
        '''
    }

}

if ( params.bowtie2 ){

    process run_bowtie2{
        tag "$pair_id"
        publishDir params.outdir, mode: 'link', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from fastpfiles_bowtie2

        output:
        set pair_id, file ('bowtie2*') into bowtie2files

        shell:
        bowtie2_threads = idv_cpu - 1
        """
        bowtie2 \
        -p ${bowtie2_threads} \
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


//start calling circRNA
if ( params.circexplorer2 && params.star ){

    process run_circexplorer2{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from starfiles
        file annotationfile
        file genomefile

        output:
        set pair_id, file ('*known.txt') into circexplorer2files

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

    process run_modify_circexplorer2{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from circexplorer2files

        output:
        set pair_id, file ('*candidates.bed') into modify_circexplorer2

        shell :
        '''
        grep circ !{query_file} \
        | grep -v chrM \
	    | awk '{print $1 "\t" $2 "\t" $3 "\t" "circexplorer2" "\t" $13 "\t" $6}' \
        > !{pair_id}_modify_circexplorer2.candidates.bed
        '''
    }

}

if ( params.ciri && params.bwa ){

    process run_ciri{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from bwafiles
        file gtffile
        file genomefile
        file ciridir

        output:
        set pair_id, file ('*.txt') into cirifiles

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

    process run_modify_ciri{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from cirifiles

        output:
        set pair_id, file ('*candidates.bed') into modify_ciri

        shell :
        '''
        cat !{query_file} \
	    | sed -e '1d' \
        | grep -v chrM \
        | awk '{print $2 "\t" $3 "\t" $4 "\t" "segemehl" "\t" $5 "\t" $11}' \
        > !{pair_id}_modify_ciri.candidates.bed
        '''
    }

}

if ( params.find_circ && params.bowtie2 ){

    process run_find_circ{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from bowtie2files
        file genomefile
        file find_circdir

        output:
        set pair_id, file ('*splice_sites.bed') into find_circfiles

        conda params.conda2dir

        shell:
        bowtie2_threads = idv_cpu - 1
        """     
        python ${find_circdir}/unmapped2anchors.py ${query_file} \
        | gzip \
        > find_circ_${pair_id}_anchors.qfa.gz

        bowtie2 \
        -p ${bowtie2_threads} \
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

    process run_modify_find_circ{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        maxForks fork_number

        input:
        set pair_id, file (query_file) from find_circfiles

        output:
        set pair_id, file ('*candidates.bed') into modify_find_circfiles

        shell :
        '''
        grep CIRCULAR !{query_file} \
        | grep -v chrM \
        | grep UNAMBIGUOUS_BP \
        | grep ANCHOR_UNIQUE \
        | awk '{print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $5 "\t" $6}' \
	    | sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 \
        > !{pair_id}_modify_finc_circ.candidates.bed
        '''
    }

}

//merge the files


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    println print_cyan( workflow.success ? "Done!" : "Oops .. something went wrong" )
    // Set up the e-mail variables
    def subject = "[nf-core/cirpipe] Successful: $workflow.runName"
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
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}






/*
 * STEP 3 - Output Description HTML
 */
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




