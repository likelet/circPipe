#!/usr/bin/env nextflow

/*
========================================================================================
                                QuantCirc
========================================================================================
 * QuantCirc was developed by Dr. Qi Zhao from Sun Yat-sen University Cancer Center.
 * it implemented a workform circRNA remapping step to quantification with a pre-build count function
 * in circletools.jar, which uses a bsj remapping bamfile and genome mapped bamfile to filter the 
 * false positive circRNA reads. 
 * Homepage / Documentation
  https://github.com/likelet/circpipe
 */



/*
 * 
 *
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>: design and implement the tools.

 * requirement 
 * hisat2 
 * bedtools 

 */


def helpMessage() {
    log.info"""
    =========================================
     QuantCirc v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the QuantCirc is as follows:

    nextflow path/to/circPipe/QuantCirc.nf --reads "path/to/*{1,2}.fq.gz" --bedfile circleRNA.bed --genomefile genome.fa --hisat2_index hisat2_index 

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --hisat2_index                Hisat2 index of whole genome sequence file name with `*.ht2` suffixed. like "GRCh38" when your files are "GRCh38.*.ht2" 
      --genomefile                  genome reference file with fasta format,  such as genome.fa .
    Options:
                                    If not set, the pipeline will create the index itself.
      --singleEnd                   Specify that the reads are single ended

    """.stripIndent()
}


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

/*
========================================================================================
                         Input file parsing
========================================================================================
*/


// set outdir 
outdir =params.outdir? params.outdir : "./" 

// input data set
genomefile = file(params.genomefile) //the genomefile
if( !genomefile.exists() ) exit 1, LikeletUtils.print_red("Missing genome file: ${genomefile}")


hisat2_index = Channel.fromPath("${params.hisat2_index}*.ht2")
            .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }



// path for placing the bedfile
bedfile= file(params.bedfile) //the genomefile
if( !bedfile.exists() ) exit 1, LikeletUtils.print_red("Missing circRNA bedfile ${params.bedfile}")





// fastqfile
(Fastpfiles_recount,Fastpfiles_hisat)=Channel.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
                                                .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
                                                .into(3) 



/*
========================================================================================
                         showing the process and files
========================================================================================
*/
log.info "INFO "+LikeletUtils.print_cyan("""
========================================================================
#    ____                        _     _____  _               
#  / __ \\                      | |   / ____|(_)          
# | |  | | _   _   __ _  _ __  | |_ | |      _  _ __  ___ 
# | |  | || | | | / _` || '_ \\ | __|| |     | || '__|/ __
# | |__| || |_| || (_| || | | || |_ | |____ | || |  | (__ 
#  \\___\\_\\ \\__,_| \\_,_ ||_| |_| \\__| \\_____||_||_| 
#
========================================================================
         """)
        .stripIndent()
log.info "INFO "+LikeletUtils.print_purple("========You are running QuantCirc with the following parameters===============")
log.info "INFO "+LikeletUtils.print_purple("Checking parameters ...")

log.info "\n"
log.info "INFO "+LikeletUtils.print_yellow("==================================Input files selected==========================")
log.info "INFO "+LikeletUtils.print_yellow("Reads :                         ") + LikeletUtils.print_green(params.reads)
log.info "INFO "+LikeletUtils.print_yellow("hisatIndex file :               ") + LikeletUtils.print_green(params.hisat2_index)
log.info "INFO "+LikeletUtils.print_yellow("Genome file :                   ") + LikeletUtils.print_green(params.genomefile)
log.info "\n"
log.info "INFO "+LikeletUtils.print_yellow("=====================================Reads types================================")
log.info "INFO "+LikeletUtils.print_yellow("SingleEnd :                     ") + LikeletUtils.print_green(params.singleEnd)
log.info "\n"
log.info "INFO "+LikeletUtils.print_yellow("==================================Output files directory========================")
log.info "INFO "+LikeletUtils.print_yellow("Output directory :              ") + LikeletUtils.print_green(outdir)
log.info "\n"




    /*
    ========================================================================================
                                    after running the tools
                                    Recount for merge
    ========================================================================================
    */
    process getPsudoCircSequenceAndBuildHisatIndex {

      storeDir "${params.outdir}/BSJ_Index"
      input:
           file bedfile
           file genomefile
      output:
           file "*.ht2" into Candidate_circRNA_index
      script:
      """
        # extract bed file for obtaining seqeuence
        sh ${baseDir}/bin/ProcessBedforGettingSequence.sh ${bedfile} temp.sort.bed temp.start.bed temp.end.bed

        bedtools getfasta -name -fi ${genomefile} -s -bed temp.start.bed > temp.start.fa
        bedtools getfasta -name -fi ${genomefile} -s -bed temp.end.bed > temp.end.fa
        # circRNA <= 400 bp
        bedtools getfasta -name -fi ${genomefile} -s -bed temp.sort.bed > temp.sort.fa 

        # merge and get combined fasta formatted psudoCirc sequences
        sh ${baseDir}/bin/MergeBSJsequence.sh temp.sort.fa temp.start.fa temp.end.fa tmp_candidate.circular_BSJ_flank.fa

        hisat2-build -p ${task.cpus}  tmp_candidate.circular_BSJ_flank.fa candidate_circRNA_BSJ_flank 
        rm temp* 
        rm tmp*
      
      """
    }

    process Recount_generate_BSJ_Bamfile {
      tag "$sampleID"
      maxForks 3
      input:
            file index from Candidate_circRNA_index.collect()
            tuple val(sampleID),  file(query_file) from Fastpfiles_recount
      output:
            tuple val(sampleID),file("${sampleID}_denovo.bam") into BSJ_mapping_bamfile
            file "fileforwaiting.txt" into Wait_for_hisat2
      script:
       if(params.singleEnd){
            """
             hisat2 -p 8 -t -k 1 -x candidate_circRNA_BSJ_flank -U ${query_file} | samtools view -bS  -q 10 -  > ${sampleID}_denovo.bam 
             touch fileforwaiting.txt
            """
        }else{
            """
            hisat2 -p 8 -t -k 1 -x candidate_circRNA_BSJ_flank -1 ${query_file[0]}  -2 ${query_file[1]} | samtools view -bS -q 10 - > ${sampleID}_denovo.bam 
            touch fileforwaiting.txt
            """
        }
    }

 process Recount_generate_genome_Bamfile {
      tag "$sampleID"
      maxForks 3
      input:
            file index from hisat2_index.collect()
            tuple val(sampleID),  file(query_file) from Fastpfiles_hisat
            file filewait from Wait_for_hisat2
      output:
            tuple val(sampleID),file("${sampleID}.bam") into Genome_remapping_bamfile

      script:
      index_base = index[0].toString() - ~/.\d.ht2/
       if(params.singleEnd){
            """
             hisat2 -p 8 -t -k 1 -x ${index_base} -U ${query_file} | samtools view -bS  -q 10 -  > ${sampleID}.bam 
            """
        }else{
            """
            hisat2 -p 8 -t -k 1 -x ${index_base} -1 ${query_file[0]}  -2 ${query_file[1]} | samtools view -bS -q 10 - > ${sampleID}.bam 
            """
        }
    }


BSJ_mapping_bamfile.combine(Genome_remapping_bamfile, by : 0 ).set{RecountBamfiles}



process Recount_estimate_step{

        publishDir "${params.outdir}/QuantCirc", pattern: "*.count", mode: 'link', overwrite: true
        input:
            tuple val(sampleID), file(bsjBamfile),file(genomeBamfile) from RecountBamfiles

            

        output:
            tuple val(sampleID),file("${sampleID}.count") into Single_sample_recount


        script:
        """
        java -jar ${baseDir}/bin/circpipetools.jar -recount -bsjbam ${bsjBamfile} -allBam ${genomeBamfile} -out ${sampleID}.count
        """
}




    // test
process Recount_results_combine{

        publishDir "${outdir}/Combination_Matrix", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
            file (query_file) from Single_sample_recount.collect()
         

        output:
            file ("multitools.exp.matrix") into (Matrix_for_circos, Plot_merge, PlotMergeCor)
        
        when:
            run_multi_tools

        script:
            """
            java -jar ${baseDir}/bin/circpipetools.jar -MM -dir ./ -suffix .count -out Quant_Circle.matrix
            """
 }


workflow.onComplete {
    log.info "INFO "+LikeletUtils.print_purple("Analysis completed ! Your result are stored in ") + LikeletUtils.print_green("${outdir}/Combination_Matrix")
}
