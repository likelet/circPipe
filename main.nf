#!/usr/bin/env nextflow

/*
 * cirPipe was implemented by Dr. Qi Zhao from Sun Yat-sen University Cancer Center.
 *
 *
 *   cirPipe is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *      See the GNU General Public License for more details.
 *
 *
 */

/*
 * to be added
 *
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>: design and implement the pipeline.
 * Wei qijin
 */

// requirement:
// - fastp/fastqc／AfterQC
// - STAR/tophat2/bowtie2/hisat2/StringTie
// - samtools/sambamba
// - Cufflinks/gffcompare
// - Bedops
// - CPAT
// - PLEK
// - CNCI
// - kallisto [https://pachterlab.github.io/kallisto/starting]

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



log.info """\
         c i r P i p e   P I P E L I N E
         =============================


         reads : ${params.reads}
         starindex : ${params.starindex}
         outdir : ${params.outdir}

         """
         .stripIndent()


/*
 * the output directory
 */
outdir = file(params.outdir)

/*
 * Add input file error exceptions Here
 */
if( !outdir.exists() ) exit 1, "Missing output directory: ${outdir}"

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

    input:
    set pair_id, file(query_file) from read_pairs_fastp

    output:
    set pair_id, file ('fastp_*') into fastpfiles_star
    set pair_id, file ('fastp_*') into fastpfiles_bwa
    file ('*.html') into fastqc_for_waiting

    """
    fastp \
    -i ${query_file[0]} \
    -I ${query_file[1]} \
    -o fastp_${pair_id}_1.fq.gz \
    -O fastp_${pair_id}_2.fq.gz
    """
}

fastqc_for_waiting = fastqc_for_waiting.first()


//start mapping
if ( params.aligner == 'star' ){
    starindex = file(params.starindex) //the index directory
    if( !starindex.exists() ) exit 1, "Missing star index directory: ${starindex}"

    process run_star{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        input:
        set pair_id, file(query_file) from fastpfiles_star
        file starindex

        output:
        set pair_id, file ('*.junction') into starfiles

        """
        STAR \
    	--runThreadN 20 \
    	--chimSegmentMin 10 \
    	--genomeDir ${starindex} \
    	--readFilesCommand zcat \
    	--readFilesIn ${query_file[0]} ${query_file[1]} \
    	--outFileNamePrefix star_${pair_id}_
        """
    }
}
else if ( params.aligner == 'bwa' ){
    process run_bwa{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        input:
        set pair_id, file (query_file) from fastpfiles_bwa

        output:
        set pair_id, file ('*.sam') into bwafiles

        """
        bwa \
    	mem -t 20 -T 19 -M -R \
    	"@RG\\tID:fastp_${pair_id}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${pair_id}" \
    	/home/wqj/test/bwaindex/genome \
    	${query_file[0]} ${query_file[1]} \
	> bwa_${pair_id}.mem.sam
        """
    }
}


//start calling circRNA
if ( params.circall == 'circexplorer2' && params.aligner == 'star' ){
    annotationfile = file(params.annotationfile) //the annotationfile
    if( !annotationfile.exists() ) exit 1, "Missing annotation file: ${annotationfile}"
    genomefile = file(params.genomefile) //the genomefile
    if( !genomefile.exists() ) exit 1, "Missing genome file: ${genomefile}"

    process run_circexplorer2{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        input:
        set pair_id, file (query_file) from starfiles
        file annotationfile
        file genomefile

        output:
        set pair_id, file ('CIRCexplorer2*') into circexplorer2files

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
}
else if ( params.circall == 'ciri' && params.aligner == 'bwa' ){
    gtffile = file(params.gtffile) //the annotationfile
    if( !gtffile.exists() ) exit 1, "Missing annotation file: ${gtffile}"
    genomefile = file(params.genomefile) //the genomefile
    if( !genomefile.exists() ) exit 1, "Missing genome file: ${genomefile}"
    ciridir = file(params.ciridir)
    if( !genomefile.exists() ) exit 1, "Missing CIRI Directory: ${ciridir}"

    process run_ciri{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        input:
        set pair_id, file (query_file) from bwafiles
        file gtffile
        file genomefile
        file ciridir

        output:
        set pair_id, file ('CIRI*') into cirifiles

        """
        perl ${ciridir}/CIRI2.pl \
	    -T 10 \
	    -F ${genomefile} \
	    -A ${gtffile} \
	    -G CIRI_${pair_id}.log \
	    -I ${query_file} \
	    -O CIRI_${pair_id}.ciri \
	    > CIRI_${pair_id}_detail.log
        """
    }
}







