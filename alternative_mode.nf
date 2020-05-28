// this nf were design for analysis from bed file and fastq file 


// input data set
genomefile = file(params.genomefile) //the genomefile
if( !genomefile.exists() ) exit 1, LikeletUtils.print_red("Missing genome file: ${genomefile}")


faifile = file(params.faifile) //the genomefile
if( !faifile.exists() ) exit 1, LikeletUtils.print_red("Missing genome file index: ${faifile}")

gtffile = file(params.gtffile) //the annotationfile-gtf-format
if( !gtffile.exists() ) exit 1, LikeletUtils.print_red("Missing gtf annotation file: ${gtffile}")

// path for placing the bedfile
inputDir= "${params.outdir}/circRNA_Identification"
if(params.bedDir){
inputDir=params.bedDir
}


// check tool number 

number_of_tools=0
String toolstring = params.selectTools 
if( toolstring.indexOf("1")!=-1){
    run_circexplorer2 = true
    number_of_tools=number_of_tools+1
}else{
    run_circexplorer2 = false
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



// fastqfile
Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { Reads_for_recount }


// circexplorer2 import
if(run_circexplorer2){
    Modify_circexplorer2 = Channel.fromPath( inputDir+'/circexplorer2_*.candidates.bed' )
    process Circexplorer2_Matrix{
        publishDir "${params.outdir}/circRNA_Identification/CIRCexplorer2", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
            file (query_file) from Modify_circexplorer2.collect()

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
        sed -i 's/_modify.candidates.bed//g' circexplorer2_merge.matrix

        echo -e "circexplorer2" > Name_circexplorer2.txt
        '''

}else{
    Merge_circexplorer2=Channel.empty()
    Name_circexplorer2=Channel.empty()
}


// ciri import
if(run_ciri){
    Modify_ciri = Channel.fromPath( inputDir+'/ciri_*.candidates.bed' )

    process Ciri_Matrix{

        publishDir "${params.outdir}/circRNA_Identification/CIRI", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
        file (query_file) from Modify_ciri.collect()


        output:
        file ('ciri_merge.matrix') into (Output_ciri,Plot_ciri_cor,Plot_ciri,Merge_ciri)
        file ('Name_ciri.txt') into Name_ciri


        shell :
        '''
        # merge sample into matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o ciri -sup 5 -merge
        mv ciri_merge.bed ciri_merge.matrix
        
        sed -i 's/ciri_//g' ciri_merge.matrix
        sed -i 's/_modify.candidates.bed//g' ciri_merge.matrix

        #for what 
        echo -e "ciri" > Name_ciri.txt
        '''
    }
}else{
    Merge_ciri=Channel.empty()
    Name_ciri=Channel.empty()
}

// Mapsplice  import 
if(run_mapsplice){
    Modify_mapsplice = Channel.fromPath( inputDir+'/mapsplice_*.candidates.bed' )

    process Mapsplice_Matrix{
        publishDir "${params.outdir}/circRNA_Identification/Mapsplice", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
        file (query_file) from Modify_mapsplice.collect()

        output:
        file ('mapsplice_merge.matrix') into (Output_mapsplice,Plot_mapsplice,Plot_mapsplice_cor,Merge_mapsplice)
        file ('Name_mapsplice.txt') into Name_mapsplice
   


        shell :
        '''
        # merge sample into matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o mapsplice -sup 5 -merge
        mv mapsplice_merge.bed mapsplice_merge.matrix

        sed -i 's/mapsplice_//g' mapsplice_merge.matrix
        sed -i 's/_modify.candidates.bed//g' mapsplice_merge.matrix
        echo -e "mapsplice" > Name_mapsplice.txt
        '''
    }

}else{
    Merge_mapsplice=Channel.empty()
    Name_mapsplice=Channel.empty()
}

// segemehl import 
if(run_segemehl){
    Modify_segemehl = Channel.fromPath( inputDir+'/segemehl_*.candidates.bed' ) 

     process Segemehl_Matrix{
            publishDir "${params.outdir}/circRNA_Identification/Segemehl", mode: 'copy', pattern: "*.matrix", overwrite: true

            input:
            file (query_file) from Modify_segemehl.collect()

            output:
            file ('segemehl_merge.matrix') into (Output_segemehl,Plot_segemehl,Merge_segemehl,Plot_segemehl_cor)
            file ('Name_segemehl.txt') into Name_segemehl
    
            shell :
            '''
            # merge sample into matrix 
            java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o segemehl -sup 5 -merge
            mv segemehl_merge.bed segemehl_merge.matrix
            # modify sample names 
            sed -i 's/segemehl_//g' segemehl_merge.matrix
            sed -i 's/_modify.candidates.bed//g' segemehl_merge.matrix
            echo -e "segemehl" > Name_segemehl.txt
            '''
    }
}else{
    Merge_segemehl=Channel.empty()
    Name_segemehl=Channel.empty()
}
// findCirc import

if(run_find_circ){
    Modify_find_circfiles = Channel.fromPath( inputDir+'/findCirc_*.candidates.bed' ) 

    process Find_circ_Matrix{
        publishDir "${params.outdir}/circRNA_Identification/Find_circ", mode: 'copy', pattern: "*.matrix", overwrite: true

        input:
        file (query_file) from Modify_find_circfiles.collect()

        output:

        file ('find_circ_merge.matrix') into (Output_find_circ,Plot_find_circ,Plot_find_circ_cor,Merge_find_circ)
        file ('Name_find_circ.txt') into Name_find_circ


        shell :
        '''

        # merge sample into matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -i candidates.bed -o find_circ -sup 5 -merge
        mv find_circ_merge.bed find_circ_merge.matrix
     
        # modify sample names 
        sed -i 's/_modify_find_circ.candidates.bed//g' find_circ_merge.matrix
        echo -e "find_circ" > Name_find_circ.txt
        '''
    }
else{
    Merge_find_circ=Channel.empty()
     Name_find_circ=Channel.empty()
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


    Combine_matrix_file= Merge_find_circ.concat( Merge_circexplorer2, Merge_ciri, Merge_mapsplice, Merge_segemehl )
    Combine_name_file=Name_find_circ.concat( Name_circexplorer2, Name_ciri, Name_mapsplice, Name_segemehl )

    /*
    ========================================================================================
                                    after running the tools
                        calculate the results by different tools
    ========================================================================================
    */
    process Tools_Merge{
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
        awk -F  "\t" '{OFS="\t"}{if ($3 > $2) {name=($1"_"$2"_"$3"_"$6);print $1,$2,$3,name,$5,$6} else {name=($1"_"$3"_"$2"_"$6);print $1,$3,$2,name,$5,$6} }' temp_concatenate.txt  | awk '$3 - $2 >= 100 && $3 - $2 <=100000 ' >  concatenate.txt
        

        for file in !{query_file}
        do 
            awk '{OFS="\t"}NR>1{print  $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t1"}' $file > ${file%%merge.matrix}merge_temp.matrix
        done 
        

        # merge and get ven merge matrix 
        java -jar !{baseDir}/bin/circpipetools.jar -collapse  -dir ./ -suffix _merge_temp.matrix -out Merged_matrix_forVen.matrix -out2 tools_merge.bed 

        awk '{OFS="\t"}{$4=".";print $0}' tools_merge.bed > all_tools_merged.matrix 
        
        awk -F  "\t" '{OFS="\t"}{if ($3 > $2) {name=($1"_"$2"_"$3"_"$6);print $1,$2,$3,name,$5,$6} else {name=($1"_"$3"_"$2"_"$6);print $1,$3,$2,name,$5,$6} }' all_tools_merged.matrix  | awk '$3 - $2 >= 100 && $3 - $2 <=100000 ' >  all_tools_merge_filtered.matrix 
        
       

      
        

        
        
        '''
    }

    

    /*
    ========================================================================================
                                    after running the tools
                                    Recount for merge
    ========================================================================================
    */
    process getPsudoCircSequenceAndBuildHisatIndex {
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

    process Recount_generate_BSJ_Bamfile {
      tag "$sampleID"
      input:
            file index from Candidate_circRNA_index.collect()
            tuple val(sampleID),  file(query_file) from Fastpfiles_recount
      output:
            tuple val(sampleID),file("${sampleID}.bam") into BSJ_mapping_bamfile
      when:
            run_multi_tools
      script:
       if(params.singleEnd){
            """
             hisat2 -p ${task.cpus} -t -k 1 -x candidate_circRNA_BSJ_flank -U ${query_file} | samtools view -bS  -q 10 -  > ${sampleID}.bam 
            """
        }else{
            """
            hisat2 -p ${task.cpus} -t -k 1 -x candidate_circRNA_BSJ_flank -1 ${query_file[0]}  -2 ${query_file[1]} | samtools view -bS -q 10 - > ${sampleID}.bam 
            """
        }
    }


    if(params.singleEnd){
        process Recount_estimate_step_single{

            input:
                tuple val(sampleID), file(bsjBamfile) from BSJ_mapping_bamfile

                

            output:
                tuple val(sampleID),file("${sampleID}.count") into Single_sample_recount

            when:
                run_multi_tools
            script:
            """
            java -jar ${baseDir}/bin/circpipetools.jar -recount -bsjbam ${bsjBamfile} -out ${sampleID}.count
            """
        }

    }else{
        process Recount_estimate_step_paired{
            tag "$sampleID"

            input:
                tuple val(sampleID), file(bsjBamfile) from BSJ_mapping_bamfile

            output:
            tuple val(sampleID),file("${sampleID}.count") into Single_sample_recount

            when:
                run_multi_tools
            script:
            """
            java -jar ${baseDir}/bin/circpipetools.jar -recount -bsjbam ${bsjBamfile} -out ${sampleID}.count --paired
            
            """
            
        }

    }

        


    // test
    process Recount_results_combine{

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


    /*
    ========================================================================================
                                        after recount
                                    Differential Expression
    ========================================================================================
    */

        process Merge_DE{
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
    process Merge_Annotation{
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
        '''
        #
        annotatePeaks.pl !{bed_file}  !{genomefile} -gtf !{gtffile} > annotated.circRNA.txt

        Rscript !{baseDir}/bin/circos.R !{baseDir}/bin/R_function.R  !{params.genomebuild} !{faifile} !{query_file}
        #perl !{baseDir}/bin/try_annotate_forGTF.pl !{gtffile} !{bed_file} newtest
        Rscript !{baseDir}/bin/circRNA_feature.R !{baseDir}/bin/R_function.R  annotated.circRNA.txt newtest.anno.txt
        '''
    }

 process Venn{
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
                                after running the tools
                                     produce report
========================================================================================
*/
process Report_production{
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
