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

params.str = 'Hello world!'

params.reads = "/home/wqj/database/test/*{1,2}.fq.gz"
params.starindex = '/home/wqj/test/starindex'
//params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
//params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
//params.outdir = 'results'

log.info """\
         c i r P i p e   P I P E L I N E
         =============================


         reads : ${params.reads}
         starindex : ${params.starindex}

         """
         .stripIndent()

/*
 * the index directory
 */
starindex = file(params.starindex)

/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }



/*
  Add input file error exceptions Here
 */

//params.input_folder = '/home/wqj/database/test'
//params.output_folder = '/home/wqj/test'


process test_fastp{
    tag "$pair_id"

    input:
    set pair_id, file(query_file) from read_pairs

    output:
    set pair_id, file ('fastp_*') into fastpfiles

    """
    fastp -i ${query_file[0]} -I ${query_file[1]} -o fastp_${pair_id}_1.fq.gz -O fastp_${pair_id}_2.fq.gz
    """
}

process test_star{
    tag "$pair_id"

    input:
    set pair_id, file(query_file) from fastpfiles
    file starindex

    output:
    set pair_id, file ('star*') into starfiles

    """
    /home/wqj/tools/STAR/bin/Linux_x86_64/STAR \
	--runThreadN 20 \
	--chimSegmentMin 10 \
	--genomeDir ${starindex} \
	--readFilesCommand zcat \
	--readFilesIn ${query_file[0]} ${query_file[1]} \
	--outFileNamePrefix star
    """
}

//not finished yet.
process test_bwa{
    input:
    file query_file from fastpfiles

    output:
    file 'out_for_bwa*' into bwafiles

    """
    /home/wqj/tools/bwa/bwa mem -t 20 -T 19 \
	-M -R "@RG\tID:out_for_fastp\tPL:PGM\tLB:noLB\tSM:out_for_fastp" \
	/home/wqj/test/bwaindex/genome \
	${query_file[0]} ${query_file[1]} > /home/wqj/test/out_for_bwa.mem.sam
    """
}
