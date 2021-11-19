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
    log.info "INFO "+"""
    =========================================
     CircPipe v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow path/to/circPipe/main.nf --reads "path/to/*{1,2}.fq.gz" -profile standard,docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --designfile                  A txt file that stored experimental design information
      --comparefile                 A txt file that stored experimental compare information

    Configuration:
      --genomefile                  Path to Fasta reference (required if not set in config file)
      --gtffile
      --annotationfile              Different annotation files from GENCODE database for annotating circRNAs. 
                                    e.g. [gencode.v25.annotation.gtf]/[gencode.v25.annotation.bed]/[hg38_gencode.txt]
      --ciridir/--find_circdir
      --mapsdir                     Home folder of ciri/find_circ/mapsplice installed location
      --genomebuild                 specific genome build for circplot, default 'hg19'; Available 'GRCh38', 'hg10'

    Options:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, test

                                    If not set, the pipeline will create the index itself.
      --singleEnd                   Specify that the reads are single ended
      --selectTools                 Specify which tools should be use. 
                                    1 for circexplorer2, 2 for ciri, 3 for find_circ, 4 for mapsplice, 5 for segemehl
                                    For example, you can set 1,2,3,4,5 for running five tools in the same time.
      --skipDE                      skip differential expression analysis                  
      --outdir                      The output directory of the results
      --mRNA                        Path to the mRNA expression matrix. Only need to be set when you want to do the correlation.


    Other options:
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

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
number_of_tools=0
String toolstring = params.selectTools 
if( toolstring.indexOf("1")!=-1){
    run_circexplorer2 = true
    number_of_tools=number_of_tools+1
}else{
    run_circexplorer2 = false
    Star_multiqc=Channel.empty()
}
if( toolstring.indexOf("2")!=-1){
    run_ciri = true
    number_of_tools=number_of_tools+1
}else{
    run_ciri = false
}
if( toolstring.indexOf("3")!=-1){
    run_find_circ = true
    number_of_tools=number_of_tools+1
}else{
    run_find_circ = false
    Bowtie2_multiqc=Channel.empty()
}
if( toolstring.indexOf("4")!=-1){
    run_mapsplice = true
    number_of_tools=number_of_tools+1
}else{
    run_mapsplice = false
}
if( toolstring.indexOf("5")!=-1){
    run_segemehl = true
    number_of_tools=number_of_tools+1
}else{
    run_segemehl = false
}


// jugde wether multiple tools involved 
run_multi_tools=false
if(number_of_tools>1){
    run_multi_tools=true
}





if(params.mRNA){
    mRNAfile = file(params.mRNA) //the mRNA file
    if( !mRNAfile.exists() ) exit 1, LikeletUtils.print_red("Missing mRNA expression file: ${params.mRNAfile}")

}


/*
========================================================================================
                         the reference directory
========================================================================================
*/




if(run_circexplorer2){
    annotationfile = file(params.annotationfile) //the annotationfile
    if( !annotationfile.exists() ) exit 1, LikeletUtils.print_red("Missing annotation file: ${params.annotationfile}")
}

genomefile = file(params.genomefile) //the genomefile
if( !genomefile.exists() ) exit 1, LikeletUtils.print_red("Missing genome file: ${params.genomefile}")

if(run_mapsplice){
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


  hisat2_index = Channel.fromPath("${params.hisat2_index}*.ht2")
            .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }

/*
========================================================================================
                         checking the design and compare file
========================================================================================
*/
//design file
if(!params.skipDE){
    designfile = file(params.designfile)
    if(params.designfile) {
        if( !designfile.exists() ) exit 1, LikeletUtils.print_red("Design file not found: ${params.designfile}")
    }
    //compare file
    comparefile = file(params.comparefile)
    if(params.comparefile){

        if( !comparefile.exists() ) exit 1, LikeletUtils.print_red("Compare file not found: ${params.comparefile}")
    }
}else {
    // add the sentence avoiding errors by nextflow 
    designfile=file(params.designfile)
    comparefile=file(params.comparefile)
}




/*
========================================================================================
                         showing the process and files
========================================================================================
*/
log.info "INFO "+LikeletUtils.print_cyan("""
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
log.info "INFO "+LikeletUtils.print_purple("============You are running circpipe with the following parameters===============")
log.info "INFO "+LikeletUtils.print_purple("Checking parameters ...")
log.info "\n"
log.info "INFO "+LikeletUtils.print_yellow("=====================================Reads types================================")
log.info "INFO "+LikeletUtils.print_yellow("SingleEnd :                     ") + LikeletUtils.print_green(params.singleEnd)
log.info "\n"
log.info "INFO "+LikeletUtils.print_yellow("====================================Tools selected==============================")
log.info "INFO "+LikeletUtils.print_yellow("Circexplorer2 :                 ") + LikeletUtils.print_green(run_circexplorer2)
log.info "INFO "+LikeletUtils.print_yellow("Find_circ :                     ") + LikeletUtils.print_green(run_find_circ)
log.info "INFO "+LikeletUtils.print_yellow("Ciri :                          ") + LikeletUtils.print_green(run_ciri)
log.info "INFO "+LikeletUtils.print_yellow("Mapsplice :                     ") + LikeletUtils.print_green(run_mapsplice)
log.info "INFO "+LikeletUtils.print_yellow("Segemehl :                      ") + LikeletUtils.print_green(run_segemehl)
log.info "\n"
log.info "INFO "+LikeletUtils.print_yellow("==================================Input files selected==========================")
log.info "INFO "+LikeletUtils.print_yellow("Reads :                         ") + LikeletUtils.print_green(params.reads)
log.info "INFO "+LikeletUtils.print_yellow("Annotation file :               ") + LikeletUtils.print_green(params.annotationfile)
log.info "INFO "+LikeletUtils.print_yellow("Genome file :                   ") + LikeletUtils.print_green(params.genomefile)
log.info "INFO "+LikeletUtils.print_yellow("GTF file :                      ") + LikeletUtils.print_green(params.gtffile)
log.info "\n"
log.info "INFO "+LikeletUtils.print_yellow("==================================Output files directory========================")
log.info "INFO "+LikeletUtils.print_yellow("Output directory :              ") + LikeletUtils.print_green(params.outdir)
log.info "INFO "+LikeletUtils.print_yellow("==================================          Others      ========================")
log.info "INFO "+LikeletUtils.print_yellow("Skip Fastqc :                   ") + LikeletUtils.print_green(params.skip_fastp)
log.info "INFO "+LikeletUtils.print_yellow("Skip DE analysis :              ") + LikeletUtils.print_green(params.skipDE)
log.info "\n"


/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { Read_pairs_fastp }



log.info "INFO "+LikeletUtils.print_yellow("===================check or build the index===============================")
/*
========================================================================================
                             check or build the index
========================================================================================
*/
//check star index
if(run_circexplorer2){
    if(params.starindex){
        starindex = Channel
                .fromPath(params.starindex)
                .ifEmpty { exit 1, "STAR index not found: ${params.starindex}" }
        star_run_index = params.starindex
    }else{
       log.info  "INFO "+LikeletUtils.print_yellow("Seems that you did not provide a STAR index for circexplorer2, circPipe will built it automaticaly. And it may take hours to prepare the reference. So you can go outside and have rest before it finished . ")
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
    starindex=Channel.empty()
    
    
}

// check bowtie2 index 
if(run_find_circ){
    if(params.bowtie2index){
        (Bowtie2index,Bowtie2index_fc)= Channel
                .fromPath(params.bowtie2index+"*.bt2")
                .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bowtie2index}, and it required by find_circ"}.into(2)


        bowtie2_run_index = params.bowtie2index
    }else{
        log.info  "INFO "+LikeletUtils.print_yellow("Seems that you did not provide a Bowtie2 index for find_circ, circPipe will built it automatically. And it may take hours to prepare the reference. So you can go outside and have rest before it finished . ")
        process makeBowtie2index {
             storeDir "${params.outdir}/reference_genome"

            input:
            file genomefile


            output:
            file "*.bt2" into Bowtie2index, Bowtie2index_fc

            script:
            """
            bowtie2-build -f ${genomefile} genome
            """
            
        }
         bowtie2_run_index = 'genome'
    }
}else{
    // avoiding throw errors  by nextflow
     (Bowtie2index,Bowtie2index_fc)=Channel.empty().into(2)
     
}

// check bowtie index
if(run_mapsplice){
    if(params.bowtieindex){
        Bowtieindex = Channel
                .fromPath(params.bowtieindex+"*.ebwt")
                .ifEmpty { exit 1, "Bowtie index not found: ${params.bowtieindex}" }
        bowtie_run_index=params.bowtieindex                
    }else{
        log.info  "INFO "+LikeletUtils.print_yellow("Seems that you did not provide a Bowtie index for mapsplice, circPipe will built it automatically. And it may take hours to prepare the reference. So you can go outside and have rest before it finished . ")
        process makeBowtieindex {
            storeDir "${params.outdir}/reference_genome"

            input:
            file genomefile

            output:
            file "*.ebwt" into Bowtieindex

            script:
            """
            bowtie-build ${genomefile} genome
            """
        }
        bowtie_run_index = "genome"
    }
}else{
    // avoiding throw errors  by nextflow
    Bowtieindex=Channel.empty()
    

}

// check bwa index 
if(run_ciri||run_multi_tools){
     LikeletUtils.print_yellow("To note that we also utilized the bwa bam file for requantification, we required that all integrated analysis run the bwa mapping step")
        
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
    bwaindex=Channel.empty()
}

// check segemehl index
if(run_segemehl){
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
    Segindex=Channel.empty()
}


log.info "INFO "+LikeletUtils.print_green("==========Index pass!...==========")
log.info "INFO "+LikeletUtils.print_green("==========Start running CircPipe...==========")



/*
========================================================================================
                       first step : run the fastp (QC tool)
========================================================================================
*/

if(params.skip_fastp){
    (Fastpfiles_mapsplice,Fastpfiles_bwa,Fastpfiles_star,Fastpfiles_segemehl,Fastpfiles_bowtie2,Fastpfiles_recount,Fastpfiles_for_sailfish,Fastpfiles_hisat)=Read_pairs_fastp.into(9)
    fastp_for_waiting=Channel.empty()
    Fastp_for_multiqc=Channel.empty()
}else{
    process RUN_FASTP{
        tag "$sampleID"
        publishDir "${params.outdir}/QC", mode: 'copy', pattern: "*_fastpreport.html", overwrite: true

        input:
        tuple val(sampleID),  file(query_file) from Read_pairs_fastp

        output:
        tuple val(sampleID),  file ('unzip_fastp_*') into (Fastqfor_swhich,Fastpfiles_mapsplice,Fastpfiles_bwa,Fastpfiles_star,Fastpfiles_segemehl,Fastpfiles_bowtie2,Fastpfiles_recount,Fastpfiles_for_sailfish,Fastpfiles_hisat)
        file ('*.html') into fastp_for_waiting
        file ('*_fastp.json') into Fastp_for_multiqc

       
        
        script:
        if(params.singleEnd){
            """
            fastp \
            -i ${query_file} \
            -o unzip_fastp_${sampleID}.fq \
            -h ${sampleID}_fastpreport.html \
            -j ${sampleID}_fastp.json
            """
        }else{
            """
            fastp \
            -i ${query_file[0]} \
            -I ${query_file[1]} \
            -o unzip_fastp_${sampleID}_1.fq \
            -O unzip_fastp_${sampleID}_2.fq \
            -h ${sampleID}_fastpreport.html \
            -j ${sampleID}_fastp.json 
            """
        }


    }
}




fastp_for_waiting = fastp_for_waiting.first() //wait for finish this process first





//   $$$$$$\  $$\                                                   $$\                                      $$$$$$\  
//  $$  __$$\ \__|                                                  $$ |                                    $$  __$$\ 
//  $$ /  \__|$$\  $$$$$$\   $$$$$$$\  $$$$$$\  $$\   $$\  $$$$$$\  $$ | $$$$$$\   $$$$$$\   $$$$$$\        \__/  $$ |
//  $$ |      $$ |$$  __$$\ $$  _____|$$  __$$\ \$$\ $$  |$$  __$$\ $$ |$$  __$$\ $$  __$$\ $$  __$$\        $$$$$$  |
//  $$ |      $$ |$$ |  \__|$$ /      $$$$$$$$ | \$$$$  / $$ /  $$ |$$ |$$ /  $$ |$$ |  \__|$$$$$$$$ |      $$  ____/ 
//  $$ |  $$\ $$ |$$ |      $$ |      $$   ____| $$  $$<  $$ |  $$ |$$ |$$ |  $$ |$$ |      $$   ____|      $$ |      
//  \$$$$$$  |$$ |$$ |      \$$$$$$$\ \$$$$$$$\ $$  /\$$\ $$$$$$$  |$$ |\$$$$$$  |$$ |      \$$$$$$$\       $$$$$$$$\ 
//   \______/ \__|\__|       \_______| \_______|\__/  \__|$$  ____/ \__| \______/ \__|       \_______|      \________|
//                                                        $$ |                                                        
//                                                        $$ |                                                        
//                                                        \__|                                                        
                                                                                                                                                                                                                                                         

if(run_circexplorer2){


    process RUN_STAR{
        tag "$sampleID"
        publishDir "${params.outdir}/Alignment/STAR", mode: 'link', overwrite: true

        input:
        tuple val(sampleID),  file(query_file) from Fastpfiles_star
        file index from starindex.collect()

        output:
        tuple val(sampleID),  file ('*.junction') into starfiles
        file ('*.out') into Star_multiqc



        shell:
        if(params.skip_fastp){
                if(params.singleEnd){
                """
                STAR \
                --runThreadN ${task.cpus} \
                --chimSegmentMin 10 \
                --genomeDir ${star_run_index} \
                --readFilesIn ${query_file} \
                --readFilesCommand zcat \
                --outFileNamePrefix star_${sampleID}_ \
                --outSAMtype None
                """
            }else{
                """
                STAR \
                --runThreadN ${task.cpus} \
                --chimSegmentMin 10 \
                --genomeDir ${star_run_index} \
                --readFilesCommand zcat \
                --readFilesIn ${query_file[0]} ${query_file[1]} \
                --outFileNamePrefix star_${sampleID}_ \
                --outSAMtype None
                """
            }
        }else{
            if(params.singleEnd){
                """
                STAR \
                --runThreadN ${task.cpus} \
                --chimSegmentMin 10 \
                --genomeDir ${star_run_index} \
                --readFilesIn ${query_file} \
                --outFileNamePrefix star_${sampleID}_ \
                --outSAMtype None
                """
            }else{
                """
                STAR \
                --runThreadN ${task.cpus} \
                --chimSegmentMin 10 \
                --genomeDir ${star_run_index} \
                --readFilesIn ${query_file[0]} ${query_file[1]} \
                --outFileNamePrefix star_${sampleID}_ \
                --outSAMtype None
                """
            }   
        }
        

    }

    /*
    ========================================================================================
                            the first tool : star - circexplorer2
                                    run the circexplorer2
    ========================================================================================
    */
    process RUN_Circexplorer2{
        tag "$sampleID"
        publishDir "${params.outdir}/circRNA_Identification/CIRCexplorer2", mode: 'copy', overwrite: true

        input:
        tuple val(sampleID),  file (query_file) from starfiles
        file annotationfile
        file genomefile

        output:
        tuple val(sampleID),  file ('*known.txt') into circexplorer2files



        script:
        """
            CIRCexplorer2 \
            parse -t STAR ${query_file} \
            > CIRCexplorer2_parse_${sampleID}.log

            CIRCexplorer2 \
            annotate -r ${annotationfile} \
            -g ${genomefile} \
            -b back_spliced_junction.bed \
            -o CIRCexplorer2_${sampleID}_circularRNA_known.txt \
            > CIRCexplorer2_annotate_${sampleID}.log
            """
    }

    /*
    ========================================================================================
                            the first tool : star - circexplorer2
                                    produce the bed6 file
    ========================================================================================
    */
    process RUN_Circexplorer2_Bed{
        tag "$sampleID"
        publishDir "${params.outdir}/circRNA_Identification/CIRCexplorer2", mode: 'copy', pattern: "*candidates.bed", overwrite: true

        input:
        tuple val(sampleID),  file (query_file) from circexplorer2files
        

        output:
        file ('*candidates.bed') into Modify_circexplorer2
        val (sampleID) into Modify_circexplorer2_id



        shell :
        '''
        if [ $((`cat !{query_file} | wc -l`)) == 0 ];then
        touch !{sampleID}_Modify_circexplorer2.candidates.bed
        else
        grep circ !{query_file} \
        | grep -v chrM \
        | awk '{print $1 "\t" $2 "\t" $3 "\t" "circexplorer2" "\t" $13 "\t" $6}' \
        > !{sampleID}_Modify_circexplorer2.temp.bed
        
        cp !{sampleID}_Modify_circexplorer2.temp.bed circexplorer2_!{sampleID}.candidates.bed
        fi
        '''
    }

    /*
    ========================================================================================
                            the first tool : star - circexplorer2
                                    produce the matrix
    ========================================================================================
    */
    process RUN_Circexplorer2_Matrix{
        publishDir "${params.outdir}/circRNA_Identification/CIRCexplorer2", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
        file (query_file) from Modify_circexplorer2.collect()
        file gtffile

        output:
        file ('circexplorer2_merge.matrix') into (Output_circexplorer2,Plot_circexplorer2,Plot_circexplorer2_cor,Merge_circexplorer2)
        file ('Name_circexplorer2.txt') into Name_circexplorer2

        shell :
        '''

        # merge sample into matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o circexplorer2 -sup 5 -merge
        mv circexplorer2_merge.bed circexplorer2_merge.matrix

        # remove non samplename string from matrix header 
        sed -i 's/circexplorer2_//g' circexplorer2_merge.matrix
        sed -i 's/.candidates.bed//g' circexplorer2_merge.matrix


        echo -e "circexplorer2" > Name_circexplorer2.txt
        '''
    }


}else{
    Merge_circexplorer2=Channel.empty()
    Name_circexplorer2=Channel.empty()
}
/*
========================================================================================
                         the first tool : star - circexplorer2
                                      run the star
========================================================================================
*/


//   $$$$$$\  $$\           $$\ 
//  $$  __$$\ \__|          \__|
//  $$ /  \__|$$\  $$$$$$\  $$\ 
//  $$ |      $$ |$$  __$$\ $$ |
//  $$ |      $$ |$$ |  \__|$$ |
//  $$ |  $$\ $$ |$$ |      $$ |
//  \$$$$$$  |$$ |$$ |      $$ |
//   \______/ \__|\__|      \__|
//                              
//                              
//                              
                                                    

/*
========================================================================================
                              the second tool : bwa - ciri
                                      run the bwa
========================================================================================
*/
if(run_ciri){
        process BWA_and_CIRI{
        tag "$sampleID"
        publishDir "${params.outdir}/circRNA_Identification/CIRI", pattern: "*.txt",mode: 'copy', overwrite: true

        input:
        tuple val(sampleID),  file (query_file) from Fastpfiles_bwa
        file index from bwaindex.collect()
        file genomefile
        file gtffile


        output:
        tuple val(sampleID),  file ('*.txt') into cirifiles

        shell:
        
        if(params.singleEnd){
            """
           bwa mem -t ${task.cpus} -T 19 -M -R "@RG\\tID:fastp_${sampleID}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${sampleID}" ${bowtie_run_index} ${query_file[0]} > ${sampleID}.sam 


            CIRI2.pl \
            -T 10 \
            -F ${genomefile} \
            -A ${gtffile} \
            -G CIRI_${sampleID}.log \
            -I ${sampleID}.sam \
            -O CIRI_${sampleID}.txt \
            > CIRI_${sampleID}_detail.log

           


             rm ${sampleID}.sam

            """
        }else{
            """
    
            bwa mem -t ${task.cpus} -T 19 -M -R "@RG\\tID:fastp_${sampleID}\\tPL:PGM\\tLB:noLB\\tSM:fastp_${sampleID}" ${bowtie_run_index} ${query_file[0]} ${query_file[1]} > ${sampleID}.sam 
            CIRI2.pl \
            -T 10 \
            -F ${genomefile}  \
            -A ${gtffile} \
            -G CIRI_${sampleID}.log \
            -I ${sampleID}.sam \
            -O CIRI_${sampleID}.txt \
            > CIRI_${sampleID}_detail.log

             # bwa sam to bam file 
           
            rm ${sampleID}.sam
            """
        }

    }


  
    /*
    ========================================================================================
                                the second tool : bwa - ciri
                                    produce the bed6 file
    ========================================================================================
    */
    process CIRI_BED{
        tag "$sampleID"
        publishDir "${params.outdir}/circRNA_Identification/CIRI", mode: 'copy', pattern: "*candidates.bed", overwrite: true

        input:
        tuple val(sampleID),  file (query_file) from cirifiles
        

        output:
        file ('*candidates.bed') into Modify_ciri

        when:
        run_ciri
        shell :
        '''
        # conver to bed file and be cautious that convert 1-based coordinate to 0-based 
            if [ $((`cat !{query_file} | wc -l`)) == 1 ];then
            touch !{sampleID}_modify_ciri.candidates.bed
            else
            cat !{query_file} \
            | sed -e '1d' \
            | grep -v chrM \
            | awk '{print $2 "\t" $3-1 "\t" $4 "\t" "ciri" "\t" $5 "\t" $11}' \
            > !{sampleID}_modify_ciri.temp.bed
            
            cp !{sampleID}_modify_ciri.temp.bed ciri_!{sampleID}.candidates.bed
            fi
            '''
    }

    /*
    ========================================================================================
                                the second tool : bwa - ciri
                                    produce the matrix
    ========================================================================================
    */
    process CIRI_Matrix{

        publishDir "${params.outdir}/circRNA_Identification/CIRI", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
        file (query_file) from Modify_ciri.collect()
        file gtffile


        output:
        file ('ciri_merge.matrix') into (Output_ciri,Plot_ciri_cor,Plot_ciri,Merge_ciri)
        file ('Name_ciri.txt') into Name_ciri


        shell :
        '''
        # merge sample into matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o ciri -sup 5 -merge
        mv ciri_merge.bed ciri_merge.matrix
        # annotate circRNA with GTFs
        # java -jar !{baseDir}/bin/circpipetools.jar -i ciri_merge.matrix -o annoted_ -gtf !{gtffile} -uniq

        sed -i 's/ciri_//g' ciri_merge.matrix
        sed -i 's/.candidates.bed//g' ciri_merge.matrix

        #for what 
        echo -e "ciri" > Name_ciri.txt
        '''
    }

    

}else{
    Merge_ciri=Channel.empty()
    Name_ciri=Channel.empty()
}



//  $$\      $$\                                         $$\ $$\                     
//  $$$\    $$$ |                                        $$ |\__|                    
//  $$$$\  $$$$ | $$$$$$\   $$$$$$\   $$$$$$$\  $$$$$$\  $$ |$$\  $$$$$$$\  $$$$$$\  
//  $$\$$\$$ $$ | \____$$\ $$  __$$\ $$  _____|$$  __$$\ $$ |$$ |$$  _____|$$  __$$\ 
//  $$ \$$$  $$ | $$$$$$$ |$$ /  $$ |\$$$$$$\  $$ /  $$ |$$ |$$ |$$ /      $$$$$$$$ |
//  $$ |\$  /$$ |$$  __$$ |$$ |  $$ | \____$$\ $$ |  $$ |$$ |$$ |$$ |      $$   ____|
//  $$ | \_/ $$ |\$$$$$$$ |$$$$$$$  |$$$$$$$  |$$$$$$$  |$$ |$$ |\$$$$$$$\ \$$$$$$$\ 
//  \__|     \__| \_______|$$  ____/ \_______/ $$  ____/ \__|\__| \_______| \_______|
//                         $$ |                $$ |                                  
//                         $$ |                $$ |                                  
//                         \__|                \__|                                  
                                                                                                                    
if(run_mapsplice){
        /*
    ========================================================================================
                                the third tool : mapsplice
                                    run the mapsplice
    ========================================================================================
    */

    process RUN_Mapsplice{
        tag "$sampleID"
        publishDir "${params.outdir}/circRNA_Identification/Mapsplice", mode: 'copy', overwrite: true

        input:

        tuple val(sampleID),  file (query_file) from Fastpfiles_mapsplice
        //tuple val(sampleID),  file (query_file) from Fastqraw_mapslice
        file gtffile
        file refmapsplice_dir from Refmapsplice.collect()
        file index from Bowtieindex.collect()

        output:
        tuple val(sampleID),  file('*') into Mapsplicefiles



        shell:
        if(params.singleEnd){
            """
            source activate find_circ
            gunzip -c ${query_file} > ${sampleID}.fastq

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
            -1 ${sampleID}.fastq \
	        --bam \
            -o output_mapsplice_${sampleID} 

            rm *.fastq
            """
        }else{
            """
            source activate find_circ
            gunzip -c ${query_file[0]} > ${sampleID}_1.fastq
            gunzip -c ${query_file[1]} > ${sampleID}_2.fastq

            mapsplice.py \
            -p ${task.cpus} \
            -k 1 \
            --fusion-non-canonical \
            --non-canonical-double-anchor \
            --min-fusion-distance 200 \
            -x ${bowtie_run_index} \
            --gene-gtf ${gtffile} \
            -c ${refmapsplice_dir} \
            -1 ${sampleID}_1.fastq \
            -2 ${sampleID}_2.fastq \
	        --bam \
            -o output_mapsplice_${sampleID} 

            rm *.fastq
            """
        }

    }

    /*
    ========================================================================================
                                the third tool : mapsplice
                                    produce the bed6 file
    ========================================================================================
    */
    process RUN_Mapsplice_Bed{
        tag "$sampleID"
        publishDir "${params.outdir}/circRNA_Identification/Mapsplice", mode: 'copy', pattern: "*candidates.bed", overwrite: true

        input:
        tuple val(sampleID),  file (query_file) from Mapsplicefiles
        file outdir
        

        output:
        file ('*candidates.bed') into Modify_mapsplice



        shell :
        '''
        if [ $((`cat output_mapsplice_!{sampleID}/circular_RNAs.txt | wc -l`)) == 0 ];then
        touch !{sampleID}_modify_mapsplice.candidates.bed
        else
        cat output_mapsplice_!{sampleID}/circular_RNAs.txt \
        | awk '{print $6}' \
        | sed -e 's/.//' \
        > !{sampleID}_mapsplice_temp1.bed

        cat output_mapsplice_!{sampleID}/circular_RNAs.txt \
        | awk '{print $1}' \
        | awk -F"~" '{print $2}' \
        > !{sampleID}_mapsplice_temp.bed

        paste !{sampleID}_mapsplice_temp.bed !{sampleID}_mapsplice_temp1.bed output_mapsplice_!{sampleID}/circular_RNAs.txt \
        | grep -v chrM \
        | awk '{if($2=="-") print $1 "\t" $4-1 "\t" $5 "\t" "mapsplice" "\t" $7 "\t" $2 ; else print $1 "\t" $5-1 "\t" $4 "\t" "mapsplice" "\t" $7 "\t" $2 }' \
        > !{sampleID}_modify_mapsplice.temp.bed
        
        cp !{sampleID}_modify_mapsplice.temp.bed mapsplice_!{sampleID}.candidates.bed
        fi
        '''
    }

    /*
    ========================================================================================
                                the third tool : mapsplice
                                    produce the matrix
    ========================================================================================
    */
    process RUN_Mapsplice_Matrix{
        publishDir "${params.outdir}/circRNA_Identification/Mapsplice", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
        file (query_file) from Modify_mapsplice.collect()
        file gtffile

        output:
        file ('mapsplice_merge.matrix') into (Output_mapsplice,Plot_mapsplice,Plot_mapsplice_cor,Merge_mapsplice)
        file ('Name_mapsplice.txt') into Name_mapsplice
   


        shell :
        '''
        # merge sample into matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o mapsplice -sup 5 -merge
        mv mapsplice_merge.bed mapsplice_merge.matrix
       
        sed -i 's/mapsplice_//g' mapsplice_merge.matrix
        sed -i 's/.candidates.bed//g' mapsplice_merge.matrix
        echo -e "mapsplice" > Name_mapsplice.txt
        '''
    }

    
   
}else{
    Merge_mapsplice=Channel.empty()
    Name_mapsplice=Channel.empty()
}





//   $$$$$$\                                                        $$\       $$\ 
//  $$  __$$\                                                       $$ |      $$ |
//  $$ /  \__| $$$$$$\   $$$$$$\   $$$$$$\  $$$$$$\$$$$\   $$$$$$\  $$$$$$$\  $$ |
//  \$$$$$$\  $$  __$$\ $$  __$$\ $$  __$$\ $$  _$$  _$$\ $$  __$$\ $$  __$$\ $$ |
//   \____$$\ $$$$$$$$ |$$ /  $$ |$$$$$$$$ |$$ / $$ / $$ |$$$$$$$$ |$$ |  $$ |$$ |
//  $$\   $$ |$$   ____|$$ |  $$ |$$   ____|$$ | $$ | $$ |$$   ____|$$ |  $$ |$$ |
//  \$$$$$$  |\$$$$$$$\ \$$$$$$$ |\$$$$$$$\ $$ | $$ | $$ |\$$$$$$$\ $$ |  $$ |$$ |
//   \______/  \_______| \____$$ | \_______|\__| \__| \__| \_______|\__|  \__|\__|
//                      $$\   $$ |                                                
//                      \$$$$$$  |                                                
//                       \______/                                                 
                                                                                                        
if(run_segemehl){
    /*
    ========================================================================================
                                the fourth tool : segemehl
                                    run the segemehl
    ========================================================================================
    */

    process RUN_Segemehl{
                tag "$sampleID"
                publishDir "${params.outdir}/circRNA_Identification/Segemehl", mode: 'copy', overwrite: true

                input:
                tuple val(sampleID),  file (query_file) from Fastpfiles_segemehl
                file genomefile
                file index from Segindex.collect()

                output:
                tuple val(sampleID),  file ("segemehl_${sampleID}.circ.sum.bed") into segemehlfiles

                shell:
                if(params.singleEnd){
                    """
                        # step 1
                    segemehl.x \
                    -t ${task.cpus} \
                    -d ${genomefile} \
                    -i ${index} \
                    -q ${query_file[0]} \
                    -S -o segemehl_${sampleID}.sam

                    # step 2
                    grep ';C;' segemehl_${sampleID}.sngl.bed  > segemehl_${sampleID}_circ_reads.bed

                    # step 3

                    haarz.x split -m 2 -q 1 \
                    -f segemehl_${sampleID}_circ_reads.bed > segemehl_${sampleID}.circ.sum.bed

                    # remove sam file for reducing storage requirement 

                    rm segemehl_${sampleID}.sam

                    """
                }else{
                    """
                    # step 1
                    segemehl.x \
                    -t ${task.cpus} \
                    -d ${genomefile} \
                    -i ${index} \
                    -q ${query_file[0]} -p ${query_file[1]}  \
                    -S -o segemehl_${sampleID}.sam

                    # step 2
                    grep ';C;' segemehl_${sampleID}.sngl.bed  > segemehl_${sampleID}_circ_reads.bed

                    # step 3

                    haarz.x split -m 2 -q 1\
                    -f segemehl_${sampleID}_circ_reads.bed > segemehl_${sampleID}.circ.sum.bed

                    # remove sam file for reducing storage requirement 

                    rm segemehl_${sampleID}.sam

                    """
                }

            }

    /*
    ========================================================================================
                                    the fourth tool : segemehl
                                        produce the bed6 file
        ========================================================================================
    */
    process RUN_Segemehl_Bed{
            tag "$sampleID"
            publishDir "${params.outdir}/circRNA_Identification/Segemehl", mode: 'copy', pattern:"*candidates.bed", overwrite: true

            input:
            set sampleID , file ( query_file ) from segemehlfiles
            

            output:
            file ('*candidates.bed') into Modify_segemehl

            shell :
            '''
            
            awk 'NR > 1 {OFS="\t";print $1,$2,$3,"segemehl",$4,$6}' segemehl_!{sampleID}.circ.sum.bed> segemehl_!{sampleID}.candidates.bed
            '''
    }

    /*
    ========================================================================================
                                    the fourth tool : segemehl
                                        produce the matrix
    ========================================================================================
    */
    process RUN_Segemehl_Matrix{
            publishDir "${params.outdir}/circRNA_Identification/Segemehl", mode: 'copy', pattern: "*.matrix", overwrite: true

            input:
            file (query_file) from Modify_segemehl.collect()
            file gtffile

            output:
            file ('segemehl_merge.matrix') into (Output_segemehl,Plot_segemehl,Merge_segemehl,Plot_segemehl_cor)
            file ('Name_segemehl.txt') into Name_segemehl
    
            shell :
            '''
            # merge sample into matrix 
            java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o segemehl -sup 5 -merge
            mv segemehl_merge.bed segemehl_merge.matrix
            sed -i 's/segemehl_//g' segemehl_merge.matrix
            sed -i 's/.candidates.bed//g' segemehl_merge.matrix
            echo -e "segemehl" > Name_segemehl.txt
            '''
    }
  

    

}else{
    Merge_segemehl=Channel.empty()
    Name_segemehl=Channel.empty()
}




//  $$$$$$$$\ $$\                 $$\                 $$\                     
//  $$  _____|\__|                $$ |                \__|                    
//  $$ |      $$\ $$$$$$$\   $$$$$$$ |       $$$$$$$\ $$\  $$$$$$\   $$$$$$$\ 
//  $$$$$\    $$ |$$  __$$\ $$  __$$ |      $$  _____|$$ |$$  __$$\ $$  _____|
//  $$  __|   $$ |$$ |  $$ |$$ /  $$ |      $$ /      $$ |$$ |  \__|$$ /      
//  $$ |      $$ |$$ |  $$ |$$ |  $$ |      $$ |      $$ |$$ |      $$ |      
//  $$ |      $$ |$$ |  $$ |\$$$$$$$ |      \$$$$$$$\ $$ |$$ |      \$$$$$$$\ 
//  \__|      \__|\__|  \__| \_______|$$$$$$\\_______|\__|\__|       \_______|
//                                    \______|                                
//                                                                            
//                                                                             

if(run_find_circ){
    process RUN_BOWTIE2{
        tag "$sampleID"
        publishDir "${params.outdir}/Alignment/Bowtie2", pattern: "*.log", mode: 'link', overwrite: true

        input:
        tuple val(sampleID),file (query_file) from Fastpfiles_bowtie2
        file index from Bowtie2index.collect()

        output:
        tuple val(sampleID), file ('bowtie2_unmapped_*') into Bowtie2files
        tuple val(sampleID), file ('bowtie2_unmapped_*') into Bowtie2files_for_autocirc
        file ('*.log') into Bowtie2_multiqc


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
            -U ${query_file} 2> bowtie2_${sampleID}.log \
            | samtools view -hbuS - \
            | samtools sort - > bowtie2_output_${sampleID}.bam

            samtools \
            view -hf 4 bowtie2_output_${sampleID}.bam \
            | samtools view -Sb - \
            > bowtie2_unmapped_${sampleID}.bam
            rm  bowtie2_output_${sampleID}.bam
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
            -2 ${query_file[1]} 2> bowtie2_${sampleID}.log \
            | samtools view -hbuS - \
            | samtools sort - > bowtie2_output_${sampleID}.bam

            samtools \
            view -hf 4 bowtie2_output_${sampleID}.bam \
            | samtools view -Sb - \
            > bowtie2_unmapped_${sampleID}.bam

            rm  bowtie2_output_${sampleID}.bam
            """
        }

    }


    /*
    ========================================================================================
                            the fifth tool : bowtie2 - find_circ
                                    run the find_circ
    ========================================================================================
    */
    process RUN_Find_circ{
        tag "$sampleID"
        publishDir "${params.outdir}/circRNA_Identification/Find_circ", mode: 'copy', overwrite: true

        input:
        tuple val(sampleID), file (query_file) from Bowtie2files
        file genomefile
        file index from Bowtie2index_fc.collect()

        output:
        tuple val(sampleID),file ('*splice_sites.bed') into find_circfiles



        shell:
        """     
         source activate find_circ
        unmapped2anchors.py ${query_file} \
        | gzip \
        > find_circ_${sampleID}_anchors.qfa.gz
        
        bowtie2 \
            -p ${task.cpus} \
            --reorder \
            --mm \
            --score-min=C,-15,0 \
            -q \
            -x ${bowtie2_run_index} \
            -U find_circ_${sampleID}_anchors.qfa.gz \
        | find_circ.py \
            -G ${genomefile} \
            -p ${sampleID}_ \
            -s find_circ_${sampleID}_stats.sites.log \
            -n find_circ \
            -R find_circ_${sampleID}_spliced_reads.fa \
            > find_circ_${sampleID}_splice_sites.bed   

        # remove temp file for reduce the usage of the disk .
        rm find_circ_${sampleID}_spliced_reads.fa 
        rm find_circ_${sampleID}_anchors.qfa.gz
        """
    }

    /*
    ========================================================================================
                            the fifth tool : bowtie2 - find_circ
                                    produce the bed6 file
    ========================================================================================
    */
    process RUN_Find_circ_Bed{
        tag "$sampleID"
        publishDir "${params.outdir}/circRNA_Identification/Find_circ", mode: 'copy', pattern: "*candidates.bed", overwrite: true

        input:
        tuple val(sampleID),  file (query_file) from find_circfiles
        

        output:
        file ('*candidates.bed') into Modify_find_circfiles

        shell :
        '''
        if [ $((`cat !{query_file} | wc -l`)) == 1 ];then
        touch !{sampleID}_modify_find_circ.candidates.bed
        else
        grep CIRCULAR !{query_file} \
        | grep -v chrM \
        | grep UNAMBIGUOUS_BP \
        | grep ANCHOR_UNIQUE \
        | awk '{print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $5 "\t" $6}' \
        | sort -t $'\t' -k 1,1 -k 2n,2 -k 3n,3 \
        > !{sampleID}_modify_find_circ.temp.bed
        
        cp !{sampleID}_modify_find_circ.temp.bed findCirc_!{sampleID}.candidates.bed
        fi
        '''
    }

    /*
    ========================================================================================
                            the fifth tool : bowtie2 - find_circ
                                    produce the matrix
    ========================================================================================
    */
    process RUN_Find_circ_Matrix{
        publishDir "${params.outdir}/circRNA_Identification/Find_circ", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
        file (query_file) from Modify_find_circfiles.collect()
        file gtffile

        output:

        file ('find_circ_merge.matrix') into (Output_find_circ,Plot_find_circ,Plot_find_circ_cor,Merge_find_circ)
        file ('Name_find_circ.txt') into Name_find_circ


        shell :
        '''

        # merge sample into matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o find_circ -sup 5 -merge
        mv find_circ_merge.bed find_circ_merge.matrix
     
        sed -i 's/findCirc_//g' find_circ_merge.matrix
        sed -i 's/.candidates.bed//g' find_circ_merge.matrix

        echo -e "find_circ" > Name_find_circ.txt
        '''
    }

 
}else{
    Merge_find_circ=Channel.empty()
     Name_find_circ=Channel.empty()
}
/*
========================================================================================
                          the fifth tool : bowtie2 - find_circ
                                   run the bowtie2
========================================================================================
*/


/*
========================================================================================
                    run the multiqc (merge the results of fastp and star)
========================================================================================
*/
process RUN_MULTIQC{
    publishDir "${params.outdir}/MultiQC", mode: 'copy', pattern: "*.html", overwrite: true

    input:
    file (query_file) from Fastp_for_multiqc.concat( Star_multiqc, Bowtie2_multiqc ).collect()

    output:
    file ('*.html') into Multiqc_results

    script:
    """
    multiqc .
    """
}



//                                     __        __                     
//                                    /  |      /  |                    
//    _______   ______   _____  ____  $$ |____  $$/  _______    ______  
//   /       | /      \ /     \/    \ $$      \ /  |/       \  /      \ 
//  /$$$$$$$/ /$$$$$$  |$$$$$$ $$$$  |$$$$$$$  |$$ |$$$$$$$  |/$$$$$$  |
//  $$ |      $$ |  $$ |$$ | $$ | $$ |$$ |  $$ |$$ |$$ |  $$ |$$    $$ |
//  $$ \_____ $$ \__$$ |$$ | $$ | $$ |$$ |__$$ |$$ |$$ |  $$ |$$$$$$$$/ 
//  $$       |$$    $$/ $$ | $$ | $$ |$$    $$/ $$ |$$ |  $$ |$$       |
//   $$$$$$$/  $$$$$$/  $$/  $$/  $$/ $$$$$$$/  $$/ $$/   $$/  $$$$$$$/ 
//                                                                      
//                                                                      
//                                                                      

if(number_of_tools==1){
    log.info LikeletUtils.print_cyan("Only one tool were selected for analysis, skip Combination  analysis Section")
    End_merge=Channel.empty()
    CorPlotMerge=Channel.empty()
    Tools_merge_html=Channel.empty()

}else{
    Combine_matrix_file= Merge_find_circ.concat( Merge_circexplorer2, Merge_ciri, Merge_mapsplice, Merge_segemehl )
    Combine_name_file=Name_find_circ.concat( Name_circexplorer2, Name_ciri, Name_mapsplice, Name_segemehl )

    /*
    ========================================================================================
                                    after running the tools
                        calculate the results by different tools
    ========================================================================================
    */
    process TOOL_MERGE{
        publishDir "${params.outdir}/Combination_Matrix", mode: 'copy', pattern: "*.matrix", overwrite: true
        
        input:
        file (query_file) from Combine_matrix_file.collect()
        file (name_file) from Combine_name_file.collect()
        file gtffile
        
        

        output:
        file ('all_tools_merge_filtered.matrix') into (Tools_merge_html,Bed_to_sailfish_cir,Bed_for_recount)
        file ('tools_merge.bed') into (Bed_for_annotation,De_merge,Cor_merge)
        file ('Merged_matrix_forVen.matrix') into Merged_file_for_Venn
        


        shell :
        '''
    
        cat *_merge.matrix >> temp_concatenate.txt

        # filtered the circRNA length less than 100bp   
        awk -F  "\t" '{OFS="\t"}{if ($3 > $2) {name=($1"_"$2"_"$3"_"$6);print $1,$2,$3,name,$5,$6} else {name=($1"_"$3"_"$2"_"$6);print $1,$3,$2,name,$5,$6} }' temp_concatenate.txt  | awk '$3 - $2 >= 100' >  concatenate.txt
        

        for file in !{query_file}
        do 
            awk '{OFS="\t"}NR>1{print  $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t1"}' $file > ${file%%merge.matrix}merge_temp.matrix
        done 
        

        # merge and get ven merge matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -collapse  -dir ./ -suffix _merge_temp.matrix -out Merged_matrix_forVen.matrix -out2 tools_merge.bed 

        awk '{OFS="\t"}{$4=".";print $0}' tools_merge.bed > all_tools_merged.matrix 
        
        awk -F  "\t" '{OFS="\t"}{if ($3 > $2) {name=($1"_"$2"_"$3"_"$6);print $1,$2,$3,name,$5,$6} else {name=($1"_"$3"_"$2"_"$6);print $1,$3,$2,name,$5,$6} }' all_tools_merged.matrix  | awk '$3 - $2 >= 100 ' >  all_tools_merge_filtered.matrix 

        
        '''
    }

    
if(!params.skipQuant){
   process getPsudoCircSequence_And_BuildHisatIndex {
      input:
           file (bed_file) from Bed_for_recount
           file genomefile
           file faifile 
      output:
           file "*.ht2" into Candidate_circRNA_index
      script:
      """
      # extract bed file for obtaining seqeuence
      sh ${baseDir}/bin/ProcessBedforGettingSequence.sh ${bed_file} temp.sort.bed temp.start.bed temp.end.bed

      bedtools getfasta -name -fi ${genomefile} -s -bed temp.start.bed > temp.start.fa
      bedtools getfasta -name -fi ${genomefile} -s -bed temp.end.bed > temp.end.fa
      # circRNA <= 400 bp
      bedtools getfasta -name -fi ${genomefile} -s -bed temp.sort.bed > temp.sort.fa 

      # merge and get combined fasta formatted psudoCirc sequences
      sh ${baseDir}/bin/MergeBSJsequence.sh temp.sort.fa temp.start.fa temp.end.fa tmp_candidate.circular_BSJ_flank.fa

      hisat2-build -p ${task.cpus}  tmp_candidate.circular_BSJ_flank.fa candidate_circRNA_BSJ_flank 
      
      """
    }

    process RECOUNT_generate_BSJ_Bamfile {
      tag "$sampleID"
      input:
            file index from Candidate_circRNA_index.collect()
            tuple val(sampleID),  file(query_file) from Fastpfiles_recount
      output:
            tuple val(sampleID),file("${sampleID}_denovo.bam") into BSJ_mapping_bamfile
            file "fileforwaiting.txt" into Wait_for_hisat2
      when:
            run_multi_tools
      script:
       if(params.singleEnd){
            """
             hisat2 -p ${task.cpus} -t -k 1 -x candidate_circRNA_BSJ_flank -U ${query_file} | samtools view -bS  -q 10 -  > ${sampleID}_denovo.bam 
             touch fileforwaiting.txt
            """
        }else{
            """
            hisat2 -p ${task.cpus} -t -k 1 -x candidate_circRNA_BSJ_flank -1 ${query_file[0]}  -2 ${query_file[1]} | samtools view -bS -q 10 - > ${sampleID}_denovo.bam 
            touch fileforwaiting.txt
            """
        }
    }

 process RECOUNT_generate_genome_Bamfile {
      tag "$sampleID"
      input:
            file index from hisat2_index.collect()
            tuple val(sampleID),  file(query_file) from Fastpfiles_hisat
            file filewait from Wait_for_hisat2
      output:
            tuple val(sampleID),file("${sampleID}.bam") into Genome_remapping_bamfile, Genome_remapping_bamfile_for_mRNAcounting
      when:
            run_multi_tools
      script:
      index_base = index[0].toString() - ~/.\d.ht2/
       if(params.singleEnd){
            """
             hisat2 -p ${task.cpus} -t -k 1 -x ${index_base} -U ${query_file} --dta | samtools view -bS  -q 10 -  > ${sampleID}.bam 
            """
        }else{
            """
            hisat2 -p ${task.cpus} -t -k 1 -x ${index_base} -1 ${query_file[0]}  -2 ${query_file[1]} --dta | samtools view -bS -q 10 - > ${sampleID}.bam 
            """
        }
    }


BSJ_mapping_bamfile.combine(Genome_remapping_bamfile, by : 0 ).set{RecountBamfiles}



if(params.singleEnd){
    process RECOUNT_estimate_step_single{

        input:
            tuple val(sampleID), file(bsjBamfile),file(genomeBamfile) from RecountBamfiles

            

        output:
            tuple val(sampleID),file("${sampleID}.count") into Single_sample_recount

        when:
            run_multi_tools
        script:
        """
        java -jar ${baseDir}/bin/circpipetools.jar -recount -bsjbam ${bsjBamfile} -allBam ${genomeBamfile} -out ${sampleID}.count
        """
    }

}else{
    process RECOUNT_estimate_step_paired{
        tag "$sampleID"

        input:
              tuple val(sampleID), file(bsjBamfile),file(genomeBamfile) from RecountBamfiles

        output:
           tuple val(sampleID),file("${sampleID}.count") into Single_sample_recount

        when:
            run_multi_tools
        script:
        """
         java -jar ${baseDir}/bin/circpipetools.jar -recount -bsjbam ${bsjBamfile} -allBam ${genomeBamfile} -out ${sampleID}.count
        
        """
        
    }

}

    


    // test
    process RECOUNT_COMBINE{

        publishDir "${params.outdir}/Combination_Matrix", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
            file (query_file) from Single_sample_recount.collect()
         

        output:
            file ("multitools.exp.matrix") into (Matrix_for_circos, Plot_merge, PlotMergeCor)
        
        when:
            run_multi_tools

        script:
            """
            java -jar ${baseDir}/bin/circpipetools.jar -MM -dir ./ -suffix .count -out multitools.exp.matrix
            """
    }

}
    /*
    ========================================================================================
                                    after running the tools
                                    Recount for merge
    ========================================================================================
    */
 

    /*
    ========================================================================================
                                        after recount
                                    Differential Expression
    ========================================================================================
    */
if(!params.skipDE){
        process MERGE_DiffernetialExpression{
        publishDir "${params.outdir}/DE_Analysis/Merge", mode: 'copy', pattern: "*", overwrite: true

        input:
        file anno_file from De_merge
        file designfile
        file comparefile
        file matrixFile from Plot_merge
        

        output:
        file ('*') into End_merge

        when:
        run_multi_tools 

        shell:
        '''
        mkdir plotdir
        Rscript !{baseDir}/bin/circRNA_DE_analysis_with_edgeR.R !{baseDir}/bin/R_function.R !{matrixFile} !{designfile} !{comparefile} !{anno_file} plotdir
        '''


    

    }
} else{
    End_merge=Channel.empty()
}

    /*
    ========================================================================================
                                        after recount
                                            Correlation
    ========================================================================================
    */

    if(params.mRNA){
        process MERGE_CorrelationAnalysis{
            publishDir "${params.outdir}/Corrrelation_Analysis/Merge", mode: 'copy', pattern: "*", overwrite: true

            input:
            file (matrix_file) from PlotMergeCor
            file (anno_file) from Cor_merge
            file mRNAfile
            
            

            when:
            run_multi_tools

            output:
            file("*") into CorPlotMerge

            shell:
            '''
            
            Rscript !{baseDir}/bin/correlation.R !{baseDir}/bin/R_function.R !{mRNAfile} !{matrix_file} !{anno_file} 
            '''
        }
    }else{
        CorPlotMerge=Channel.empty()
    }


    /*
    ========================================================================================
                                    after running the tools
                                            annotation
    ========================================================================================
    */
    process ANNOTATION{
        publishDir "${params.outdir}/Annotation", mode: 'copy', pattern: "*", overwrite: true

        input:
        file (bed_file) from Bed_for_annotation
        file (query_file) from Matrix_for_circos
        file gtffile 
        file faifile
        file genomefile

        when:
        run_multi_tools

        output:
        file ('*') into Annotation_plot

        shell:
        """
        java -jar ${baseDir}/bin/circpipeTools.jar -annote -bedfile !{bed_file}  -gtf !{gtffile} -out annotated.circRNA.txt

        Rscript !{baseDir}/bin/circos.R !{baseDir}/bin/R_function.R  !{params.genomebuild} !{faifile} !{query_file}
        Rscript !{baseDir}/bin/circRNA_feature.R !{baseDir}/bin/R_function.R  annotated.circRNA.txt newtest.anno.txt

        """
        
    }

 process VENN_PLOT{
        publishDir "${params.outdir}/Annotation", mode: 'copy', pattern: "*", overwrite: true

        input:
        file (matrix_file) from Merged_file_for_Venn
        


        when:
        (number_of_tools < 6) && (number_of_tools > 1)

        output:
        file ('*') into venn_plot

        shell:
        '''
        Rscript !{baseDir}/bin/venn.R !{matrix_file} venn.png
        '''

}
   
}


/*
========================================================================================
                                mRNA abundence counting
========================================================================================
*/

if(params.mRNA){
    process mRNAmeasurementByFeatureCount {
      input:
        file gtffile
        tuple val(sampleID),file(bamfile)  from  Genome_remapping_bamfile_for_mRNAcounting
      output:
        tuple val(sampleID),file("${sampleID}.gene.count") into Feature_count_out
      script:
      """
      
        featureCounts -a $gtf -T ${task.cpus} -p -t gene -g gene_id -o ${sampleID}.gene.count ${bamfile}

      """
    }
    
    process featureCount2Matrix {
      publishDir "${params.outdir}/mRNAexpression", mode: 'copy', pattern: "*", overwrite: true
      input:
        file bamfiles  from  Genome_remapping_bamfile_for_mRNAcounting.collect()
      output:
        file "featureCount.matrix" into Expression_matrix_of_mRNA
      script:
      """
        java -jar ${baseDir}/bin/circpipeTools.jar -MF -dir ./ -out featureCount.matrix
      """
    }

}




/*
========================================================================================
                                after running the tools
                                     produce report
========================================================================================
*/
process REPORT_GENERATION{
    publishDir "${params.outdir}/Report", mode: 'copy', pattern: "*.html", overwrite: true

    input:
    file (de_file) from End_merge.collect()
    file (cor_file) from CorPlotMerge.collect()
    file (anno_file) from Annotation_plot.collect()
    file (calculate_file) from Tools_merge_html
    file (multiqc_file) from Multiqc_results
    
    

    when:
    run_multi_tools

    output:
    file ('*.html') into report_html

    shell:
    '''
    ln -s !{baseDir}/bin/*.Rmd ./
    Rscript -e "require( 'rmarkdown' ); render('report.Rmd', 'html_document')"
    '''
}


/*
* Completion e-mail notification
*/
if(params.email){
    emailaddress = params.email
}




workflow.onComplete {

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
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Circexplorer2 : ${run_circexplorer2} </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:cneter; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Find_circ : ${run_find_circ} </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Ciri : ${run_ciri} </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Mapsplice : ${run_mapsplice} </pre></td>
    </tr>
    <tr>
    <td style = 'text-align:center; padding: 8px; line-height: 1.42857143; vertical-align: top; border-top: 1px solid #ddd;' ><pre style="white-space: pre-wrap; overflow: visible;"> Segemehl : ${run_segemehl} </pre></td>
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






