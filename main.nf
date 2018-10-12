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
 * Add input file error exceptions Here
 */
outdir = file(params.outdir) //the output directory
if( !outdir.exists() ) exit 1, "Missing output directory: ${outdir}"

if ( params.aligner == 'star' ){

    starindex = file(params.starindex) //the index directory
    if( !starindex.exists() ) exit 1, "Missing star index directory: ${starindex}"

}
else if ( params.aligner == 'mapsplice' ){

    condadir = file(params.condadir) //the python3 environment
    if( !condadir.exists() ) exit 1, "Missing python3 environment: ${condadir}"

    gtffile = file(params.gtffile) //the annotationfile
    if( !gtffile.exists() ) exit 1, "Missing annotation file: ${gtffile}"

    mapsdir = file(params.mapsdir) //the mapsplice directory
    if( !mapsdir.exists() ) exit 1, "Missing Mapsplice Directory: ${mapsdir}"

    refdir = file(params.refdir) //the reference genome directory
    if( !refdir.exists() ) exit 1, "Missing Reference Genome Directory: ${refdir}"

}
else if ( params.aligner == 'segemehl' ){

    condadir = file(params.condadir) //the python3 environment
    if( !condadir.exists() ) exit 1, "Missing python3 environment: ${condadir}"

    genomefile = file(params.genomefile) //the genomefile
    if( !genomefile.exists() ) exit 1, "Missing genome file: ${genomefile}"

    segdir = file(params.mapsdir) //the segemehl directory
    if( !segdir.exists() ) exit 1, "Missing Segemehl Directory: ${segdir}"

    segindex = file(params.segindex) //the segemehl index file
    if( !segindex.exists() ) exit 1, "Missing Segemehl index file: ${segindex}"

}

if ( params.circall == 'circexplorer2' ){

    annotationfile = file(params.annotationfile) //the annotationfile
    if( !annotationfile.exists() ) exit 1, "Missing annotation file: ${annotationfile}"

    genomefile = file(params.genomefile) //the genomefile
    if( !genomefile.exists() ) exit 1, "Missing genome file: ${genomefile}"

}
else if ( params.circall == 'ciri' ){

    gtffile = file(params.gtffile) //the annotationfile
    if( !gtffile.exists() ) exit 1, "Missing annotation file: ${gtffile}"

    genomefile = file(params.genomefile) //the genomefile
    if( !genomefile.exists() ) exit 1, "Missing genome file: ${genomefile}"

    ciridir = file(params.ciridir)
    if( !genomefile.exists() ) exit 1, "Missing CIRI Directory: ${ciridir}"

}
else if ( params.circall == 'find_circ' ){

    conda2dir = file(params.conda2dir) //the python2 environment
    if( !conda2dir.exists() ) exit 1, "Missing python2 environment: ${conda2dir}"

    genomefile = file(params.genomefile) //the genomefile
    if( !genomefile.exists() ) exit 1, "Missing genome file: ${genomefile}"

    find_circdir = file(params.find_circdir)
    if( !find_circdir.exists() ) exit 1, "Missing find_circ Directory: ${find_circdir}"

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

    input:
    set pair_id, file(query_file) from read_pairs_fastp

    output:
    set pair_id, file ('fastp_*') into fastpfiles_star
    set pair_id, file ('fastp_*') into fastpfiles_bwa
    set pair_id, file ('fastp_*') into fastpfiles_mapsplice
    set pair_id, file ('fastp_*') into fastpfiles_segemehl
    set pair_id, file ('fastp_*') into fastpfiles_bowtie2
    file ('*.html') into fastqc_for_waiting

	shell:
	if ( params.aligner == 'mapsplice' || params.aligner == 'segemehl' ){
	    """
        fastp \
        -i ${query_file[0]} \
        -I ${query_file[1]} \
        -o fastp_${pair_id}_1.fq \
        -O fastp_${pair_id}_2.fq
        """
	} else {
        """
        fastp \
        -i ${query_file[0]} \
        -I ${query_file[1]} \
        -o fastp_${pair_id}_1.fq.gz \
        -O fastp_${pair_id}_2.fq.gz
        """
	}

}

fastqc_for_waiting = fastqc_for_waiting.first()


//start mapping
if ( params.aligner == 'star' ){

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
else if ( params.aligner == 'mapsplice' ){

    process run_mapsplice{
        tag "$pair_id"
	    publishDir params.outdir, mode: 'copy', overwrite: true

	    input:
	    set pair_id, file (query_file) from fastpfiles_mapsplice
	    file mapsdir
	    file gtffile
	    file refdir
        file outdir

	    output:
        set pair_id, file ('*.log') into mapsplicefiles

        conda params.condadir

	    """
	    python ${mapsdir}/mapsplice.py \
	    -p 25 -k 1 \
	    --fusion-non-canonical \
	    --non-canonical-double-anchor \
	    --min-fusion-distance 200 \
	    -x /home/wqj/test/bowtieindex/chrX \
	    --gene-gtf ${gtffile} \
	    -c ${refdir} \
	    -1 ${query_file[0]} \
	    -2 ${query_file[1]} \
	    -o ${outdir}/output_mapsplice_${pair_id} 2\
	    > mapsplice_${pair_id}.log
	    """
    }

}
else if ( params.aligner == 'segemehl' ){

    process run_segemehl{
        tag "$pair_id"
	    publishDir params.outdir, mode: 'copy', overwrite: true

	    input:
	    set pair_id, file (query_file) from fastpfiles_segemehl
	    file segdir
	    file genomefile
	    file segindex

	    output:
        set pair_id, file ('segemehl*') into segemehlfiles

        conda params.condadir

	    """
	    ${segdir}/segemehl.x \
	    -d ${genomefile} \
	    -i ${segindex} \
	    -q ${query_file[0]} \
	    -p ${query_file[1]} \
	    -t 20 -S \
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
else if ( params.aligner == 'bowtie2' ){

    process run_bowtie2{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        input:
        set pair_id, file (query_file) from fastpfiles_bowtie2

        output:
        set pair_id, file ('bowtie2*') into bowtie2files

        """
        bowtie2 \
        -p 20 \
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
if ( params.circall == 'circexplorer2' && params.aligner == 'star' ){

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
else if ( params.circall == 'find_circ' && params.aligner == 'bowtie2' ){

    process run_find_circ{
        tag "$pair_id"
        publishDir params.outdir, mode: 'copy', overwrite: true

        input:
        set pair_id, file (query_file) from bowtie2files
        file genomefile
        file find_circdir

        output:
        set pair_id, file ('find_circ*') into find_circfiles

        conda params.conda2dir

        """
        python ${find_circdir}/unmapped2anchors.py ${query_file} \
        | gzip \
        > find_circ_${pair_id}_anchors.qfa.gz

        bowtie2 \
        -p 20 \
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
        -n ${pair_id} \
        -R find_circ_${pair_id}_spliced_reads.fa \
        > find_circ_${pair_id}_splice_sites.bed

        grep CIRCULAR find_circ_${pair_id}_splice_sites.bed \
        | grep -v chrM \
        | awk '\$5>=2' \
        | grep UNAMBIGUOUS_BP \
        | grep ANCHOR_UNIQUE \
        | python ${find_circdir}/maxlength.py 100000 \
        > finc_circ_${pair_id}.candidates.bed
        """
    }
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}



