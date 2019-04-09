#!/usr/bin/env nextflow

/*
========================================================================================
                            circPipe
========================================================================================
 * circPipe was implemented by Dr. Qi Zhao and Qijin Wei from Sun Yat-sen University Cancer Center.
 * Homepage / Documentation
  https://github.com/likelet/circpipe

 */



/*
 * to be added
 *
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>: design and implement the pipeline.
 * Wei qijin
 */


def helpMessage() {
    log.info"""
    =========================================
     CircPipe v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow path/to/circPipe/main.nf --reads="path/to/*{1,2}.fq.gz" -profile standard,docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --designfile                  A txt file that stored experimental design information
      --comparefile                 A txt file that stored experimental compare information

    Configuration:
      --genomefile                  Path to Fasta reference (required if not set in config file)
      --gtffile/
      --annotationfile              Different annotation files from GENCODE database for annotating circRNAs. 
                                    e.g. [gencode.v25.annotation.gtf]/[gencode.v25.annotation.bed]/[hg38_gencode.txt]
      --ciridir/--find_circdir/
      --mapsdir/--knifedir          Home folder of ciri/find_circ/mapsplice/knife installed location

    Options:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

                                    If not set, the pipeline will create the index itself.
      --singleEnd                   Specify that the reads are single ended
      --merge                       Merge the different matrixes produce by different tools and draw the venn graph
      --separate                    Annotate the results separately
      --selectTools                 Specify which tools should be use. 
                                    1 for circexplorer2, 2 for ciri, 3 for find_circ, 4 for mapsplice, 5 for segemehl, 6 for knife. 
                                    For example, you can set 1,2,3,4,5 for running five tools in the same time.
      --outdir                      The output directory of the results
      --mRNA                        Path to the mRNA expression matrix. Only need to be set when you want to do the correlation.


    Other options:
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
 * Checking the input files
 * Build up the output files
 * Adding input files error exceptions Here
 */


outdir = file(params.outdir) //the output directory


// check variables 
assert params.library == 'rnase' || params.aligner == 'total' : "Invalid library option: ${params.library}. Valid options: 'rnase', 'total'"


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







if(params.mRNA){
    mRNAfile = file(params.mRNA) //the mRNA file
    if( !mRNAfile.exists() ) exit 1, LikeletUtils.print_red("Missing mRNA expression file: ${mRNAfile}")

}


/*
========================================================================================
                         the reference directory
========================================================================================
*/




if(params.circexplorer2){
    annotationfile = file(params.annotationfile) //the annotationfile
    if( !annotationfile.exists() ) exit 1, LikeletUtils.print_red("Missing annotation file: ${annotationfile}")
}

genomefile = file(params.genomefile) //the genomefile
if( !genomefile.exists() ) exit 1, LikeletUtils.print_red("Missing genome file: ${genomefile}")

if(params.mapsplice){
    if(params.refmapsplice){
        Refmapsplice = Channel.fromPath(params.refmapsplice) //the mapsplice reference genome directory
    }else{
        process built_refmapsplice_reference_by_split {
                storeDir "${params.outdir}/reference_genome/split"
                input:
                file genomefile
                output: 
                file "split" into Refmapsplice
                shell:
                """
                mkdir split 
                perl ${baseDir}/bin/split_fasta_by_chromsome.pl ${genomefile} split
                """
            }
    }
    
}else{
     Refmapsplice = Channel.fromPath(params.inputdir)
}


faifile = file(params.faifile) //the genomefile
if( !faifile.exists() ) exit 1, LikeletUtils.print_red("Missing genome file index: ${faifile}")

gtffile = file(params.gtffile) //the annotationfile-gtf-format
if( !gtffile.exists() ) exit 1, LikeletUtils.print_red("Missing gtf annotation file: ${gtffile}")




/*
========================================================================================
                         checking the design and compare file
========================================================================================
*/
//design file
designfile = file(params.designfile)
if(params.designfile) {
    if( !designfile.exists() ) exit 1, LikeletUtils.print_red("Design file not found: ${params.designfile}")
}
//compare file
comparefile = file(params.comparefile)
if(params.comparefile){

    if( !comparefile.exists() ) exit 1, LikeletUtils.print_red("Compare file not found: ${params.comparefile}")
}



/*
========================================================================================
                         showing the process and files
========================================================================================
*/
log.info LikeletUtils.print_cyan("""
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
log.info LikeletUtils.print_purple("============You are running circpipe with the following parameters===============")
log.info LikeletUtils.print_purple("Checking parameters ...")
log.info "\n"
log.info LikeletUtils.print_yellow("=====================================Reads types================================")
log.info LikeletUtils.print_yellow("SingleEnd :                     ") + LikeletUtils.print_green(params.singleEnd)
log.info "\n"
log.info LikeletUtils.print_yellow("====================================Tools selected==============================")
log.info LikeletUtils.print_yellow("Circexplorer2 :                 ") + LikeletUtils.print_green(params.circexplorer2)
log.info LikeletUtils.print_yellow("Find_circ :                     ") + LikeletUtils.print_green(params.find_circ)
log.info LikeletUtils.print_yellow("Ciri :                          ") + LikeletUtils.print_green(params.ciri)
log.info LikeletUtils.print_yellow("Mapsplice :                     ") + LikeletUtils.print_green(params.mapsplice)
log.info LikeletUtils.print_yellow("Segemehl :                      ") + LikeletUtils.print_green(params.segemehl)
log.info LikeletUtils.print_yellow("Knife :                         ") + LikeletUtils.print_green(params.knife)
log.info "\n"
log.info LikeletUtils.print_yellow("==================================Input files selected==========================")
log.info LikeletUtils.print_yellow("Reads :                         ") + LikeletUtils.print_green(params.reads)
log.info LikeletUtils.print_yellow("Annotation file :               ") + LikeletUtils.print_green(params.annotationfile)
log.info LikeletUtils.print_yellow("Genome file :                   ") + LikeletUtils.print_green(params.genomefile)
log.info LikeletUtils.print_yellow("GTF file :                      ") + LikeletUtils.print_green(params.gtffile)
log.info "\n"
log.info LikeletUtils.print_yellow("==================================Output files directory========================")
log.info LikeletUtils.print_yellow("Output directory :              ") + LikeletUtils.print_green(params.outdir)
log.info LikeletUtils.print_yellow("==================================          Others      ========================")
log.info LikeletUtils.print_yellow("Skip Fastqc :                   ") + LikeletUtils.print_green(params.skip_fastp)
log.info "\n"


/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { read_pairs_fastp }



log.info LikeletUtils.print_yellow("===================check or build the index===============================")
/*
========================================================================================
                             check or build the index
========================================================================================
*/
if(params.circexplorer2){
    if(params.starindex){
        starindex = Channel
                .fromPath(params.starindex)
                .ifEmpty { exit 1, "STAR index not found: ${params.starindex}" }
        star_run_index = params.starindex
    }else{
        LikeletUtils.print_yellow("Seems that you did not provide a STAR index for circexplorer2, circPipe will built it automaticaly. And it may take hours to prepare the reference. So you can go outside and have rest before it finished . ")
        process makeSTARindex {
            storeDir "${params.outdir}/reference_genome"

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
        star_run_index = 'starindex'
    }
}else{
    // avoiding throw errors  by nextflow
    starindex=Channel.create()
    
}


if(params.find_circ){
    if(params.bowtie2index){
        (Bowtie2index,Bowtie2index_fc)= Channel
                .fromPath(params.bowtie2index+"*.bt2")
                .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bowtie2index}, and it required by find_circ"}.into(2)


        bowtie2_run_index = params.bowtie2index
    }else{
        LikeletUtils.print_yellow("Seems that you did not provide a Bowtie2 index for find_circ, circPipe will built it automatically. And it may take hours to prepare the reference. So you can go outside and have rest before it finished . ")
        process makeBowtie2index {
             storeDir "${params.outdir}/reference_genome"

            input:
            file genomefile


            output:
            file "*.bt2" into Bowtie2index, Bowtie2index_fc,Bowtie2_build_knife

            script:
            """
            bowtie2-build -f ${genomefile} genome
            """
            
        }
         bowtie2_run_index = 'genome'
    }
}else{
    // avoiding throw errors  by nextflow
     (Bowtie2index,Bowtie2index_fc)=Channel.create().into(2)
        
}

if(params.mapsplice){
    if(params.bowtieindex){
        Bowtieindex = Channel
                .fromPath(params.bowtieindex+"*.ebwt")
                .ifEmpty { exit 1, "Bowtie index not found: ${params.bowtieindex}" }
        bowtie_run_index=params.bowtieindex                
    }else{
         LikeletUtils.print_yellow("Seems that you did not provide a Bowtie index for mapsplice, circPipe will built it automatically. And it may take hours to prepare the reference. So you can go outside and have rest before it finished . ")
        process makeBowtieindex {
            storeDir "${params.outdir}/reference_genome"

            input:
            file genomefile

            output:
            file "*.ebwt" into Bowtieindex,Bowtie_build_knife

            script:
            """
            bowtie-build ${genomefile} genome
            """
        }
        bowtie_run_index = "genome"
    }
}else{
    // avoiding throw errors  by nextflow
    (Bowtieindex,Bowtie_build_knife)=Channel.create().into(2)
}


if(params.ciri){
    if(params.bwaindex){
        bwaindex = Channel
                .fromPath(params.bwaindex+"*.{ann,amb,pac,bwt,sa}")
                .ifEmpty { exit 1, "BWA index not found: ${params.bwaindex}" }
        bowtie_run_index = params.bwaindex
    }else{
        LikeletUtils.print_yellow("Seems that you did not provide a BWA index for ciri, circPipe will built it automatically. And it may take hours to prepare the reference. So you can go outside and have rest before it finished . ")
        process makeBWAindex {
            storeDir "${params.outdir}/reference_genome"

            input:
            file genomefile


            output:
            file "*.{ann,amb,pac,bwt,sa}" into bwaindex

            script:
            """
             bwa index -p genome ${genomefile} 
            """
            bowtie_run_index="genome"
        }

    }
}else{
    // avoiding throw errors  by nextflow
    bwaindex=Channel.create()
}

if(params.segemehl){
    if(params.segindex){
       Segindex = Channel
                .fromPath(params.segindex+"*.idx")
                .ifEmpty { exit 1, "Segemehl index not found: ${params.segindex}" }
    }else{
        LikeletUtils.print_yellow("Seems that you did not provide a segemehl index for runing segemehl, circPipe will built it automatically. And it may take hours to prepare the reference. So you can go outside and have rest before it finished . ")
        process makeSegemehlindex {
            storeDir "${params.outdir}/reference_genome"

            input:
            file genomefile
            file segdir

            output:
            file "genome.idx" into Segindex

            script:
            """
        segemehl.x -d ${genomefile} -x genome.idx
        """
        }
    }
}else{
    // avoiding throw errors  by nextflow
    Segindex=Channel.create()
}




log.info LikeletUtils.print_green("==========Index pass!...==========")
log.info LikeletUtils.print_green("==========Start running CircPipe...==========")



/*
========================================================================================
                       first step : run the fastp (QC tool)
========================================================================================
*/

if(params.skip_fastp){
    (fastpfiles_mapsplice,fastpfiles_bwa,fastpfiles_star,fastpfiles_segemehl,fastpfiles_knife,Fastpfiles_bowtie2,Fastpfiles_recount)=read_pairs_fastp.into(7)
    fastp_for_waiting=Channel.create()
    fastp_for_multiqc=Channel.create()
}else{
process Fastp{
    tag "$pair_id"
    publishDir "${params.outdir}/QC", mode: 'copy', pattern: "*_fastpreport.html", overwrite: true

    input:
    set pair_id, file(query_file) from read_pairs_fastp

    output:
    set pair_id, file ('unzip_fastp_*') into fastpfiles_mapsplice,fastpfiles_bwa,fastpfiles_star,fastpfiles_segemehl,fastpfiles_knife,Fastpfiles_bowtie2,Fastpfiles_recount
    file ('*.html') into fastp_for_waiting
    file ('*_fastp.json') into fastp_for_multiqc

    when:
    !params.skip_fastp
    
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
}




fastp_for_waiting = fastp_for_waiting.first() //wait for finish this process first




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
    file ('*.out') into star_multiqc

    when:
    params.circexplorer2

    shell:
    if(params.singleEnd){
        """
        STAR \
        --runThreadN ${task.cpus} \
        --chimSegmentMin 10 \
        --genomeDir ${star_run_index} \
        --readFilesIn ${query_file} \
        --outFileNamePrefix star_${pair_id}_
        """
    }else{
        """
        STAR \
        --runThreadN ${task.cpus} \
        --chimSegmentMin 10 \
        --genomeDir ${star_run_index} \
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
    

    output:
    file ('*candidates.bed') into modify_circexplorer2
    val (pair_id) into modify_circexplorer2_id

    when:
    params.circexplorer2

    shell :
    '''
    if [ $((`cat !{query_file} | wc -l`)) == 0 ];then
    touch !{pair_id}_modify_circexplorer2.candidates.bed
    else
    grep circ !{query_file} \
    | grep -v chrM \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" "circexplorer2" "\t" $13 "\t" $6}' \
    > !{pair_id}_modify_circexplorer2.temp.bed
    
    python !{baseDir}/bin/quchong.py !{pair_id}_modify_circexplorer2.temp.bed circexplorer2_!{pair_id}_modify.candidates.bed
    fi
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
    
    file designfile
    file gtffile

    output:
    file ('circexplorer2.txt') into merge_circexplorer2
    file ('*.matrix') into output_circexplorer2
    file ('name_circexplorer2.txt') into name_circexplorer2
    file ('*annote.txt') into de_circexplorer2
    file ('*.matrix') into plot_circexplorer2
    file ('*annote.txt') into cor_circexplorer2
    file ('*.matrix') into plot_circexplorer2_cor

    when:
    params.circexplorer2

    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> concatenate.bed
    done
    
    python !{baseDir}/bin/hebinglist.py concatenate.bed merge_concatenate.bed
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
    cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
    cat mergeconcatenate.bed > circexplorer2.txt
        
    cat circexplorer2.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{baseDir}/bin/bed1114.jar -i annotation.bed -o circexplorer2_ -gtf !{gtffile} -uniq

    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_circexplorer2.candidates.bed | wc -l`)) == 0 ];then
        python !{baseDir}/bin/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{baseDir}/bin/quchongsamples.py mergeconcatenate.bed ${line}_modify_circexplorer2.candidates.bed counts.txt
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
    
    file designfile
    file comparefile
    file (matrix_file) from plot_circexplorer2
    

    output:
    file ('*') into end_circexplorer2

    when:
    params.circexplorer2 && params.separate

    shell:
    '''
    Rscript !{baseDir}/bin/edgeR_circ.R !{baseDir}/bin/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}

/*
========================================================================================
                           the first tool : star - circexplorer2
                                        Correlation
========================================================================================
*/
if(params.mRNA){
    process Circexplorer2_Cor{
        publishDir "${params.outdir}/Corrrelation_Analysis/Circexplorer2", mode: 'copy', pattern: "*", overwrite: true

        input:
        file (matrix_file) from plot_circexplorer2_cor
        file (anno_file) from cor_circexplorer2
        file mRNAfile
        
        

        when:
        params.mRNA && params.circexplorer2 && params.separate

        output:
        file ('*') into cor_plot_circexplorer2

        shell:
        '''
        Rscript !{baseDir}/bin/correlation.R !{baseDir}/bin/R_function.R !{mRNAfile} !{matrix_file} !{anno_file}
        '''
    }
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
    file genomefile


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
        ${bowtie_run_index} \
        ${query_file} \
        > bwa_${pair_id}.mem.sam
        """
    }else{
        """
        bwa \
        mem -t ${task.cpus} \
        -T 19 -M -R \
        "@RG\\tID:fastp_${pair_id}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${pair_id}" \
        ${bowtie_run_index} \
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

    output:
    set pair_id, file ('*.txt') into cirifiles

    when:
    params.ciri

    script:
        """
        CIRI2.pl \
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
    

    output:
    file ('*candidates.bed') into modify_ciri_file
    val (pair_id) into modify_ciri_id

    when:
    params.ciri

    shell :
    '''
        if [ $((`cat !{query_file} | wc -l`)) == 1 ];then
        touch !{pair_id}_modify_ciri.candidates.bed
        else
        cat !{query_file} \
	    | sed -e '1d' \
        | grep -v chrM \
        | awk '{print $2 "\t" $3 "\t" $4 "\t" "ciri" "\t" $5 "\t" $11}' \
        > !{pair_id}_modify_ciri.temp.bed
        
        python !{baseDir}/bin/quchong.py !{pair_id}_modify_ciri.temp.bed !{pair_id}_modify_ciri.candidates.bed
        fi
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
    
    file designfile
    file gtffile


    output:
    file ('ciri.txt') into merge_ciri
    file ('*.matrix') into output_ciri
    file ('name_ciri.txt') into name_ciri
    file ('*annote.txt') into de_ciri
    file ('*.matrix') into plot_ciri
    file ('*annote.txt') into cor_ciri
    file ('*.matrix') into plot_ciri_cor

    when:
    params.ciri

    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> concatenate.bed
    done
    
    python !{baseDir}/bin/hebinglist.py concatenate.bed merge_concatenate.bed
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
    cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
    cat mergeconcatenate.bed > ciri.txt    
        
    cat ciri.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{baseDir}/bin/bed1114.jar -i annotation.bed -o ciri_ -gtf !{gtffile} -uniq

    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_ciri.candidates.bed | wc -l`)) == 0 ];then
        python !{baseDir}/bin/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{baseDir}/bin/quchongsamples.py mergeconcatenate.bed ${line}_modify_ciri.candidates.bed counts.txt
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
    
    file designfile
    file comparefile
    file (matrix_file) from plot_ciri
    

    output:
    file ('*') into end_ciri

    when:
    params.ciri && params.separate

    shell:
    '''
    Rscript !{baseDir}/bin/edgeR_circ.R !{baseDir}/bin/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}

/*
========================================================================================
                                the second tool : bwa - ciri
                                        Correlation
========================================================================================
*/
if(params.mRNA){
    process Ciri_Cor{
        publishDir "${params.outdir}/Corrrelation_Analysis/CIRI", mode: 'copy', pattern: "*", overwrite: true

        input:
        file (matrix_file) from plot_ciri_cor
        file (anno_file) from cor_ciri
        file mRNA
        
        

        when:
        params.mRNA && params.ciri && params.separate

        output:
        file ('*') into cor_plot_ciri

        shell:
        '''
        Rscript !{baseDir}/bin/correlation.R !{baseDir}/bin/R_function.R !{mRNA} !{matrix_file} !{anno_file}
        '''
    }
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
    file gtffile
    file refmapsplice_dir from Refmapsplice
    file outdir
    file index from Bowtieindex.collect()

    output:
    set pair_id, file('*') into mapsplicefiles


    when:
    params.mapsplice

    shell:
    if(params.singleEnd){
        """
        mapsplice.py \
        -p ${task.cpus} \
        -k 1 \
        --fusion-non-canonical \
        --qual-scale phred33 \
        --non-canonical-double-anchor \
        --min-fusion-distance 200 \
        -x ${bowtie_run_index} \
        --gene-gtf ${gtffile} \
        -c ${refmapsplice_dir} \
        -1 ${query_file} \
        -o output_mapsplice_${pair_id} 2 \
        > ${pair_id}_mapsplice.log
        """
    }else{
        """
        mapsplice.py \
        -p ${task.cpus} \
        -k 1 \
        --fusion-non-canonical \
        --non-canonical-double-anchor \
        --min-fusion-distance 200 \
        -x ${bowtie_run_index} \
        --gene-gtf ${gtffile} \
        -c ${refmapsplice_dir} \
        -1 ${query_file[0]} \
        -2 ${query_file[1]} \
        -o output_mapsplice_${pair_id} 2 \
        > ${pair_id}_mapsplice.log
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
    
    python !{baseDir}/bin/quchong.py !{pair_id}_modify_mapsplice.temp.bed mapsplice_!{pair_id}_modify.candidates.bed
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
    
    file designfile
    file gtffile

    output:
    file ('mapsplice.txt') into merge_mapsplice
    file ('*.matrix') into output_mapsplice
    file ('name_mapsplice.txt') into name_mapsplice
    file ('*annote.txt') into de_mapsplice
    file ('*.matrix') into plot_mapsplice
    file ('*annote.txt') into cor_mapsplice
    file ('*.matrix') into plot_mapsplice_cor

    when:
    params.mapsplice

    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> concatenate.bed
    done
    
    python !{baseDir}/bin/hebinglist.py concatenate.bed merge_concatenate.bed
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
    cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
    cat mergeconcatenate.bed > mapsplice.txt
        
    cat mapsplice.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{baseDir}/bin/bed1114.jar -i annotation.bed -o mapsplice_ -gtf !{gtffile} -uniq
    
    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_mapsplice.candidates.bed | wc -l`)) == 0 ];then
        python !{baseDir}/bin/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{baseDir}/bin/quchongsamples.py mergeconcatenate.bed ${line}_modify_mapsplice.candidates.bed counts.txt
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
    
    file designfile
    file comparefile
    file (matrix_file) from plot_mapsplice
    

    output:
    file ('*') into end_mapsplice

    when:
    params.mapsplice && params.separate

    shell:
    '''
    Rscript !{baseDir}/bin/edgeR_circ.R !{baseDir}/bin/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}

/*
========================================================================================
                                  the third tool : mapsplice
                                        Correlation
========================================================================================
*/
if(params.mRNA){
    process Mapsplice_Cor{
        publishDir "${params.outdir}/Corrrelation_Analysis/Mapsplice", mode: 'copy', pattern: "*", overwrite: true

        input:
        file (matrix_file) from plot_mapsplice_cor
        file (anno_file) from cor_mapsplice
        file mRNAfile
        
        

        when:
        params.mapsplice && params.separate

        output:
        file ('*') into cor_plot_mapsplice

        shell:
        '''
        Rscript !{baseDir}/bin/correlation.R !{baseDir}/bin/R_function.R !{mRNAfile} !{matrix_file} !{anno_file}
        '''
    }
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
            file genomefile
            file index from Segindex.collect()

            output:
            set pair_id, file ('*splicesites.bed') into segemehlfiles

            when:
            params.segemehl 

            shell:
            if(params.singleEnd){
                """
                segemehl.x \
                -d ${genomefile} \
                -i ${params.segindex} \
                -q ${query_file} \
                -t ${task.cpus} \
                -S \
                | samtools view -bS - \
                | samtools sort -o - \
                | samtools view -h - \
                > segemehl_${pair_id}_mapped.sam

                testrealign.x \
                -d ${genomefile} \
                -q segemehl_${pair_id}_mapped.sam \
                -n \
                -U segemehl_${pair_id}_splicesites.bed \
                -T segemehl_${pair_id}_transrealigned.bed
                """
            }else{
                """
                segemehl.x \
                -d ${genomefile} \
                 -i ${params.segindex} \
                -q ${query_file[0]} \
                -p ${query_file[1]} \
                -t ${task.cpus} \
                -S \
                | samtools view -bS - \
                | samtools sort -o - \
                | samtools view -h - \
                > segemehl_${pair_id}_mapped.sam

                testrealign.x \
                -d ${genomefile} \
                -q segemehl_${pair_id}_mapped.sam \
                -n \
                -U segemehl_${pair_id}_splicesites.bed \
                -T segemehl_${pair_id}_transrealigned.bed
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
        

        output:
        file ('*candidates.bed') into modify_segemehl
        val (pair_id) into modify_segemehl_id

        when:
        params.segemehl 

        shell :
        '''
        if [ $((`cat !{query_file} | wc -l`)) == 0 ];then
        touch !{pair_id}_modify_segemehl.candidates.bed
        else
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
        
        python !{baseDir}/bin/quchong.py !{pair_id}_modify_segemehl.temp.bed segemehl_!{pair_id}_modify.candidates.bed
        fi
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
        
        file designfile
        file gtffile

        output:
        file ('segemehl.txt') into merge_segemehl
        file ('*.matrix') into output_segemehl
        file ('*annote.txt') into de_segemehl
        file ('*.matrix') into plot_segemehl
        file ('name_segemehl.txt') into name_segemehl
        file ('*annote.txt') into cor_segemehl
        file ('*.matrix') into plot_segemehl_cor


        when:
        params.segemehl 

        shell :
        '''
        for file in !{query_file}
        do
            cat $file >> concatenate.bed
        done
        
        python !{baseDir}/bin/hebinglist.py concatenate.bed merge_concatenate.bed
        sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
        cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
        cat mergeconcatenate.bed > segemehl.txt
        
        cat segemehl.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
        java -jar !{baseDir}/bin/bed1114.jar -i annotation.bed -o segemehl_ -gtf !{gtffile} -uniq

        cat !{designfile} > designfile.txt
        sed -i '1d' designfile.txt
        cat designfile.txt | awk '{print $1}' > samplename.txt
        
        echo -e "id\\c" > merge_header.txt
        
        cat samplename.txt | while read line
        do
            if [ $((`cat ${line}_modify_segemehl.candidates.bed | wc -l`)) == 0 ];then
            python !{baseDir}/bin/createzero.py mergeconcatenate.bed counts.txt
            else
            python !{baseDir}/bin/quchongsamples.py mergeconcatenate.bed ${line}_modify_segemehl.candidates.bed counts.txt
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
        
        file designfile
        file comparefile
        file (matrix_file) from plot_segemehl
        

        output:
        file ('*') into end_segemehl

        when:
        params.separate && params.segemehl 

        shell:
        '''
        Rscript !{baseDir}/bin/edgeR_circ.R !{baseDir}/bin/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
        '''
    }

    /*
    ========================================================================================
                                    the fourth tool : segemehl
                                            Correlation
    ========================================================================================
    */
    if(params.mRNA){
        process Segemehl_Cor{
            publishDir "${params.outdir}/Corrrelation_Analysis/Segemehl", mode: 'copy', pattern: "*", overwrite: true

            input:
            file (matrix_file) from plot_segemehl_cor
            file (anno_file) from cor_segemehl
            file mRNAfile
            
            

            when:
            params.segemehl && params.separate

            output:
            file ('*') into cor_plot_segemehl

            shell:
            '''
            Rscript !{baseDir}/bin/correlation.R !{baseDir}/bin/R_function.R !{mRNAfile} !{matrix_file} !{anno_file}
            '''
        }
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
    set pair_id, file (query_file) from Fastpfiles_bowtie2
    file index from Bowtie2index.collect()

    output:
    set pair_id, file ('bowtie2_unmapped_*') into Bowtie2files
    set pair_id, file ('bowtie2_unmapped_*') into Bowtie2files_for_autocirc
    file ('*.log') into bowtie2_multiqc

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
        -x ${bowtie2_run_index} \
        -q \
        -U ${query_file} 2> bowtie2_${pair_id}.log \
        | samtools view -hbuS - \
        | samtools sort - > bowtie2_output_${pair_id}.bam

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
        -x ${bowtie2_run_index} \
        -q \
        -1 ${query_file[0]} \
        -2 ${query_file[1]} 2> bowtie2_${pair_id}.log \
        | samtools view -hbuS - \
        | samtools sort - > bowtie2_output_${pair_id}.bam

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
    set pair_id, file (query_file) from Bowtie2files
    file genomefile
    file index from Bowtie2index_fc.collect()

    output:
    set pair_id, file ('*splice_sites.bed') into find_circfiles


    when:
    params.find_circ

    shell:
    """     
    unmapped2anchors.py ${query_file} \
    | gzip \
    > find_circ_${pair_id}_anchors.qfa.gz

    bowtie2 \
    -p ${task.cpus} \
    --reorder \
    --mm \
    --score-min=C,-15,0 \
    -q \
    -x ${bowtie2_run_index} \
    -U find_circ_${pair_id}_anchors.qfa.gz \
    | find_circ.py \
    -G ${genomefile} \
    -p ${pair_id}_ \
    -s find_circ_${pair_id}_stats.sites.log \
    -n find_circ \
    -R find_circ_${pair_id}_spliced_reads.fa \
    > find_circ_${pair_id}_splice_sites.bed   
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
    

    output:
    file ('*candidates.bed') into modify_find_circfiles
    val (pair_id) into modify_find_circ_id

    when:
    params.find_circ

    shell :
    '''
    if [ $((`cat !{query_file} | wc -l`)) == 1 ];then
    touch !{pair_id}_modify_find_circ.candidates.bed
    else
    grep CIRCULAR !{query_file} \
    | grep -v chrM \
    | grep UNAMBIGUOUS_BP \
    | grep ANCHOR_UNIQUE \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $5 "\t" $6}' \
    | sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 \
    > !{pair_id}_modify_find_circ.temp.bed
    
    python !{baseDir}/bin/quchong.py !{pair_id}_modify_find_circ.temp.bed !{pair_id}_modify_find_circ.candidates.bed
    fi
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
    
    python !{baseDir}/bin/hebinglist.py concatenate.bed merge_concatenate.bed
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
    cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
    cat mergeconcatenate.bed > find_circ.txt
    
    cat find_circ.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{baseDir}/bin/bed1114.jar -i annotation.bed -o find_circ_ -gtf !{gtffile} -uniq

    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_find_circ.candidates.bed | wc -l`)) == 0 ];then
        python !{baseDir}/bin/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{baseDir}/bin/quchongsamples.py mergeconcatenate.bed ${line}_modify_find_circ.candidates.bed counts.txt
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
    
    file designfile
    file comparefile
    file (matrix_file) from plot_find_circ
    

    output:
    file ('*') into end_find_circ

    when:
    params.find_circ && params.separate

    shell:
    '''
    Rscript !{baseDir}/bin/edgeR_circ.R !{baseDir}/bin/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}

/*
========================================================================================
                          the fifth tool : bowtie2 - find_circ
                                        Correlation
========================================================================================
*/
if(params.mRNA){
    process Find_circ_Cor{
        publishDir "${params.outdir}/Corrrelation_Analysis/Find_circ", mode: 'copy', pattern: "*", overwrite: true

        input:
        file (matrix_file) from plot_find_circ_cor
        file (anno_file) from cor_find_circ
        file mRNAfile
        
        

        when:
        params.find_circ && params.separate

        output:
        file ('*') into cor_plot_find_circ

        shell:
        '''
        Rscript !{baseDir}/bin/correlation.R !{baseDir}/bin/R_function.R !{mRNAfile} !{matrix_file} !{anno_file}
        '''
    }
}


/*
========================================================================================
                              the sixth tool : knife
                                   run the knife
========================================================================================
*/
process Knife{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/KNIFE", mode: 'copy', overwrite: true

    input:
    set pair_id, file (query_file) from fastpfiles_knife



    output:
    set pair_id, file ('*.txt') into knifefiles


    when:
    params.knife

    shell:
    knifedir = params.knifedir
    knfieprefix = params.knifeprefix
    if(params.singleEnd){
        '''
        ln -s !{knifedir} ./
        
        pwd | awk '{print "sh completeRun.sh",$0,"complete",$0,"testData 15", !{knfieprefix}, "circReads 40 1 2>&1"}' | bash

        cp ./testData/circReads/combinedReports/naiveunzip* ./
        '''
    }else{
        '''
        ln -s !{knifedir} ./
       
        pwd | awk '{print "sh completeRun.sh",$0,"appended",$0,"testData 13", !{knfieprefix}, "circReads 40 1 2>&1"}' | bash
        
        cp ./testData/circReads/combinedReports/naiveunzip* ./
        '''
    }

}

/*
========================================================================================
                              the sixth tool : knife
                              produce the bed6 file
========================================================================================
*/
process Knife_Bed{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/KNIFE", mode: 'copy', pattern:"*candidates.bed", overwrite: true

    input:
    set pair_id , file ( query_file ) from knifefiles
    

    output:
    file ('*candidates.bed') into modify_knifefiles
    val (pair_id) into modify_knife_id

    when:
    params.knife

    shell :
    '''
    if [ $((`cat !{query_file} | grep -v reg | wc -l`)) == 1 ];then
    touch !{pair_id}_modify_knife.candidates.bed
    else
    cat !{query_file} \
    | awk 'NR>1{print $10}' \
    > reads.txt

    cat !{query_file} \
    | awk 'NR>1{print $1}' \
    | awk -F"|" '{print $1}' \
    > chr.txt

    cat !{query_file} \
    | awk 'NR>1{print $1}' \
    | awk -F"|" '{print $5}' \
    > strand.txt

    cat !{query_file} \
    | awk 'NR>1{print $1}' \
    | awk -F"|" '{print $2}' \
    | awk -F":" '{print $2}' \
    > start.txt

    cat !{query_file} \
    | awk 'NR>1{print $1}' \
    | awk -F"|" '{print $3}' \
    | awk -F":" '{print $2}' \
    > end.txt

    cat !{query_file} \
    | awk 'NR>1{print $1}' \
    | awk -F"|" '{print $4}' \
    > kind.txt

    paste chr.txt start.txt end.txt reads.txt strand.txt kind.txt \
    | grep -v reg \
    | grep -v chrM \
    | awk '{if($5=="-") print $1 "\t" $2 "\t" $3 "\t" "knife" "\t" $4 "\t" $5 ; else print $1 "\t" $3 "\t" $2 "\t" "knife" "\t" $4 "\t" $5 }' \
    > temp.bed

    python !{baseDir}/bin/quchong.py temp.bed !{pair_id}_modify_knife.candidates.bed
    fi
    '''

}

/*
========================================================================================
                               the sixth tool : knife
                                 produce the matrix
========================================================================================
*/
process Knife_Matrix{
    tag "$pair_id"
    publishDir "${params.outdir}/circRNA_Identification/KNIFE", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
    file (query_file) from modify_knifefiles.collect()
    val (pair_id) from modify_knife_id.collect()
    
    file designfile
    file gtffile

    output:
    file ('knife.txt') into merge_knife
    file ('*annote.txt') into de_knife
    file ('*.matrix') into output_knife
    file ('*.matrix') into plot_knife
    file ('name_knife.txt') into name_knife
    file ('*annote.txt') into cor_knife
    file ('*.matrix') into plot_knife_cor


    when:
    params.knife

    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> concatenate.bed
    done
    
    python !{baseDir}/bin/hebinglist.py concatenate.bed merge_concatenate.bed
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 merge_concatenate.bed > mergeconcatenate.bed 
    cat mergeconcatenate.bed | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id.txt
    cat mergeconcatenate.bed > knife.txt
    
    cat knife.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > annotation.bed
    java -jar !{baseDir}/bin/bed1114.jar -i annotation.bed -o knife_ -gtf !{gtffile} -uniq

    cat !{designfile} > designfile.txt
    sed -i '1d' designfile.txt
    cat designfile.txt | awk '{print $1}' > samplename.txt
    
    echo -e "id\\c" > merge_header.txt
    
    cat samplename.txt | while read line
    do
        if [ $((`cat ${line}_modify_knife.candidates.bed | wc -l`)) == 0 ];then
        python !{baseDir}/bin/createzero.py mergeconcatenate.bed counts.txt
        else
        python !{baseDir}/bin/quchongsamples.py mergeconcatenate.bed ${line}_modify_knife.candidates.bed counts.txt
        fi
        paste -d"\t" id.txt counts.txt > temp.txt
        cat temp.txt > id.txt
        echo -e "\\t${line}\\c" >> merge_header.txt
    done   
    
    sed -i 's/\\[//g' merge_header.txt
    sed -i 's/\\,//g' merge_header.txt
    sed -i 's/\\]//g' merge_header.txt
    echo -e "\\n\\c" >> merge_header.txt
     
    cat merge_header.txt id.txt > knife_merge.matrix
    echo -e "knife" > name_knife.txt
    '''
}

/*
========================================================================================
                                 the sixth tool : knife
                                 Differential Expression
========================================================================================
*/
process Knife_DE{
    publishDir "${params.outdir}/DE_Analysis/Knife", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (anno_file) from de_knife
    
    file designfile
    file comparefile
    file (matrix_file) from plot_knife
    

    output:
    file ('*') into end_knife

    when:
    params.knife && params.separate

    shell:
    '''
    Rscript !{baseDir}/bin/edgeR_circ.R !{baseDir}/bin/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}

/*
========================================================================================
                                  the sixth tool : knife
                                        Correlation
========================================================================================
*/
if(params.mRNA){
    process Knife_Cor{
        publishDir "${params.outdir}/Corrrelation_Analysis/Knife", mode: 'copy', pattern: "*", overwrite: true

        input:
        file (matrix_file) from plot_knife_cor
        file (anno_file) from cor_knife
        file mRNAfile
        
        

        when:
         params.knife && params.separate

        output:
        file ('*') into cor_plot_knife

        shell:
        '''
        Rscript !{baseDir}/bin/correlation.R !{baseDir}/bin/R_function.R !{mRNAfile} !{matrix_file} !{anno_file}
        '''
    }
}


/*
========================================================================================
                    run the multiqc (merge the results of fastp and star)
========================================================================================
*/
process Multiqc{
    publishDir "${params.outdir}/MultiQC", mode: 'copy', pattern: "*.html", overwrite: true

    input:
    file (query_file) from fastp_for_multiqc.concat( star_multiqc, bowtie2_multiqc ).collect()

    output:
    file ('*.html') into multiqc_results

    script:
    """
    multiqc .
    """
}



/*
========================================================================================
                                after running the tools
                       calculate the results by different tools
========================================================================================
*/
process Tools_Merge{
    publishDir "${params.outdir}/Combination_Matrix", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
    file (query_file) from merge_find_circ.concat( merge_circexplorer2, merge_ciri, merge_mapsplice, merge_segemehl, merge_knife ).collect()
    file (name_file) from name_find_circ.concat( name_circexplorer2, name_ciri, name_mapsplice, name_segemehl, name_knife ).collect()
    
    

    output:
    file ('all_tools_merge.matrix') into tools_merge
    file ('all_tools_merge.matrix') into tools_merge_html
    file ('for_annotation.bed') into Med_for_annotation
    file ('for_annotation.bed') into Bed_for_recount
    file ('for_annotation.bed') into bed_for_merge
    file ('*annote.txt') into de_merge
    file ('*annote.txt') into cor_merge
    file ('all_tools_intersect.matrix') into tools_intersect


    shell :
    '''
    for file in !{query_file}
    do
        cat $file >> temp_concatenate.txt
    done 
        
    awk '$3-$2>=100' temp_concatenate.txt > concatenate.txt
    
    python !{baseDir}/bin/hebingtoolsid.py concatenate.txt id_unsort.txt
    sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 id_unsort.txt > id_sort.txt
    
    cat id_sort.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $4 }' > for_annotation.bed
    java -jar !{baseDir}/bin/bed1114.jar -i for_annotation.bed -o merge_ -gtf !{gtffile} -uniq 
    
    echo -e "total\\c" > total.txt
    cat id_sort.txt | awk '{print $1 "_" $2 "_" $3 "_" $4 }' > id_merge.txt
    cat id_merge.txt total.txt > id_list.txt
    
    for file in !{query_file}
    do
        python !{baseDir}/bin/countnumbers.py id_sort.txt $file counts.txt
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
    
    Rscript !{baseDir}/bin/intersect.R all_tools_merge.matrix
    '''
}


/*
========================================================================================
                                after running the tools
                                   Recount for merge
========================================================================================
*/
process Recount_index_step{
    input:
        file (bed_file) from Bed_for_recount
        file genomefile
        file faifile

    output:
        file "*.ht2" into Candidate_circRNA_index
        file ('tmp_candidate_circRNA.gff3') into Gff3_file

    when:
        params.merge

    shell:
        '''
        # sort bed (in some result bed file , start > end ) and length filtering( >= 100nt)
        awk -F  "\t" '{OFS="\t"}{if ($3 > $2) {name=($1"_"$2"_"$3"_"$6);print $1,$2,$3,name,$5,$6} else {name=($1"_"$3"_"$2"_"$6);print $1,$3,$2,name,$5,$6} }' !{bed_file} | awk '$3 - $2 >= 100 ' >  tmp_candidate_circRNA.bed

        # bed to gff3 for htseq-count; sites around junction sites(+/-3bp)
        awk  '{OFS="\t"}{split($4,a,"_");len=$3-$2; print $4"("a[4]")",".","exon",len-3,len+3,".","+",".","gene_id="$4 }' tmp_candidate_circRNA.bed > tmp_candidate_circRNA.gff3

        # bed to fasta
        bedtools getfasta -fi  !{genomefile} -s -bed tmp_candidate_circRNA.bed -name > tmp_candidate.circular.fa

        # candidate circRNA sequnces (doulbed).
        awk 'NR%2==1{print $0}NR%2==0{print $1$1}'  tmp_candidate.circular.fa > tmp_candidate.circular_doulbed.fa

        #build index for candidate circRNA sequnce
        mkdir $tmp_candidate_hisat_index

        hisat2-build -p !{task.cpus} tmp_candidate.circular_doulbed.fa candidate_circRNA_doulbed 



        '''
}

process Recount_estimate_step{

    input:
        file index from Candidate_circRNA_index.collect()
        file (gff_file) from Gff3_file
        set pair_id, file(query_file) from Fastpfiles_recount
        

    output:
        file('*circRNA_requantity.count') into single_sample_recount

    when:
        params.merge

    shell:
    if(params.singleEnd){
        '''
        sh !{baseDir}/bin/final_recount2.sh !{pair_id} single_end !{task.cpus} !{query_file}
        '''
    }else{
        '''
        sh !{baseDir}/bin/final_recount2.sh !{pair_id} pair_end !{task.cpus} !{query_file[0]} !{query_file[1]}
        '''
    }
}


// test
process Recount_results_combine{

    publishDir "${params.outdir}/Combination_Matrix", mode: 'copy', pattern: "*.matrix", overwrite: true

    input:
        file (query_file) from single_sample_recount.collect()
        file designfile
        file (bed_file) from bed_for_merge

    output:
        file ("final.matrix") into (Matrix_for_circos, Plot_merge, PlotMergeCor)
    
    when:
        params.merge

    shell:
        '''
        cat !{designfile} > designfile.txt
        sed -i '1d' designfile.txt
        cat designfile.txt | awk '{print $1}' > samplename.txt
        
        echo -e "id\\c" > merge_header.txt
        
        cat for_annotation.bed | awk '{print $1 "_" $2 "_" $3 "_" $6 }' > id.txt
        
        cat samplename.txt | while read line
        do
            sed '$d' ${line}_circRNA_requantity.count > temp1.bed
            sed '$d' temp1.bed > temp2.bed
            sed '$d' temp2.bed > temp3.bed
            sed '$d' temp3.bed > temp4.bed
            sed '$d' temp4.bed > ${line}_modify_circRNA_requantity.count
            python !{baseDir}/bin/final_countnumbers.py id.txt ${line}_modify_circRNA_requantity.count ${line}_counts.txt
            paste -d"\t" id.txt ${line}_counts.txt > temp.txt
            cat temp.txt > id.txt
            echo -e "\\t${line}\\c" >> merge_header.txt
        done   
        
        echo -e "\\n\\c" >> merge_header.txt
        
        cat merge_header.txt id.txt > final.matrix
        '''
}


/*
========================================================================================
                                      after recount
                                 Differential Expression
========================================================================================
*/
process Merge_DE{
    publishDir "${params.outdir}/DE_Analysis/Merge", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (anno_file) from de_merge
    
    file designfile
    file comparefile
    file (matrix_file) from Plot_merge
    

    output:
    file ('*') into end_merge

    when:
    params.merge

    shell:
    '''
    Rscript !{baseDir}/bin/edgeR_circ.R !{baseDir}/bin/R_function.R !{matrix_file} !{designfile} !{comparefile} !{anno_file}
    '''
}

/*
========================================================================================
                                       after recount
                                        Correlation
========================================================================================
*/

if(params.mRNA){
    process Merge_Cor{
        publishDir "${params.outdir}/Corrrelation_Analysis/Merge", mode: 'copy', pattern: "*", overwrite: true

        input:
        file (matrix_file) from PlotMergeCor
        file (anno_file) from cor_merge
        file mRNAfile
        
        

        when:
        params.merge

        output:
        file("*") into CorPlotMerge

        shell:
        '''
        Rscript !{baseDir}/bin/correlation.R !{baseDir}/bin/R_function.R !{mRNAfile} !{matrix_file} !{anno_file}
        '''
    }
}else{
    CorPlotMerge=Channel.create()
}


/*
========================================================================================
                                after running the tools
                                        annotation
========================================================================================
*/
process Merge_Annotation{
    publishDir "${params.outdir}/Annotation", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (bed_file) from Med_for_annotation
    file (query_file) from Matrix_for_circos
    file gtffile 

    when:
    params.merge

    output:
    file ('*') into annotation_plot

    shell:
    '''
    java -jar !{baseDir}/bin/bed1114.jar -i !{bed_file} -o merge_ -gtf !{gtffile} -uniq 
    Rscript !{baseDir}/bin/circos.R !{query_file}
    perl !{baseDir}/bin/try_annotate_forGTF.pl !{gtffile} !{bed_file} newtest
    Rscript !{baseDir}/bin/circRNA_feature.R !{baseDir}/bin/R_function.R merge_for_annotation_annote.txt newtest.anno.txt
    '''
}

process Venn{
    publishDir "${params.outdir}/Annotation", mode: 'copy', pattern: "*", overwrite: true

    input:
    file (matrix_file) from tools_merge
    


    when:
    params.selectTools=='1,2,3,4,6'

    output:
    file ('*') into venn_plot

    shell:
    '''
    Rscript !{baseDir}/bin/venn.R !{matrix_file} venn.png
    '''
}


/*
========================================================================================
                                after running the tools
                                     produce report
========================================================================================
*/
process Report_production{
    publishDir "${params.outdir}/Report", mode: 'copy', pattern: "*.html", overwrite: true

    input:
    file (de_file) from end_merge.collect()
    file (cor_file) from CorPlotMerge.collect()
    file (anno_file) from annotation_plot.collect()
    file (calculate_file) from tools_merge_html
    file (multiqc_file) from multiqc_results
    
    

    when:
    params.merge

    output:
    file ('*.html') into report_html

    shell:
    '''
    cp !{baseDir}/bin/*.Rmd ./
    Rscript -e "require( 'rmarkdown' ); render('report.Rmd', 'html_document')"
    '''
}


/*
* Completion e-mail notification
*/

emailaddress = params.email



workflow.onComplete {

println LikeletUtils.print_cyan( workflow.success ? "Done!" : "Oops .. something went wrong" )

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
    <h1> CircPipe execution summary </h1>
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
        <pre style="white-space: pre-wrap; overflow: visible; margin-bottom: 0;">${workflow.errorMessage}</pre>
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
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Knife : ${params.knife} </pre></td>
    </tr>
    </tbody>
    </table>

    <h4> likelet/CircPipe </h4>
    <h4><a href="https://github.com/likelet/circPipe">https://github.com/likelet/circPipe</a></h4>
    <h4> If you need help, you can send email to Qi Zhao(zhaoqi@sysucc.org.cn) or Wei Qijin (513848731@qq.com) </h4>
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
/*   def subject = "[circpipe] Successful: $workflow.runName"
if(!workflow.success){
    subject = "[circpipe] FAILED: $workflow.runName"
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


log.info "[circpipe] Pipeline Complete"
*/



