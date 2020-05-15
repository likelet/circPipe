# circPipe
**circRNA analysis pipe**

[![Build Status](https://travis-ci.org/likelet/cirpipe.svg?branch=master)](https://travis-ci.org/likelet/cirpipe)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![Docker](https://img.shields.io/docker/automated/likelet/cirpipe.svg)](https://hub.docker.com/r/likelet/cirpipe)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
In the recent decades, circular RNA(circRNA) have came into our sight and emerged as a large class of non-coding RNAs. It thought to be predominately generated by covalently back-splicing conjunction of 3’end and 5’end from the same or distinct exons. circRNA involving in a variaty of biological functions in cells. Emerging evidences show that dysregulation of circRNAs are associated with several human diseases including cancer. Therefore, Systematically identification and profiling circRNAs is a fundimental process for  dissecting the underlining biological mechanism of its regulation. To date, a range of tools have been developed to investigate circRNAs from high-throughput sequencing including CIRCexplorer2, CIRI, Find_circ, Mapsplice, and Segemehl. However, it not a easy thing for users to pick tools and compare the result from such tools. The most appropriate strategy is run those tools parallel and collapse the result together, and vote the candidates to improve their confidence for further experimental exploration. In addition, a remapping step is nessasary to help reduce the false postive candidates in an application scenarios. Here, we present circPipe, a nextflow-based pipeline for running multitool-based identification of circRNA from RNA-seq dataset. circPipe integrates a remapping stragegy to help filter out the circRNAs with non-reads supported in non-mismatch mode of realignment. The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation   

A full tutorial of CircPipe can be found at Wiki page of this project. plz go to the https://github.com/likelet/circPipe/wiki


### Pipeline Steps

CircPipe allows you to run arbitrary pipelines among five pipelines.
Choose between workflows by using `--selectTools` or not(default) .

| Step                                         | Pipeline One     | Pipeline Two          | Pipeline Three        | Pipeline Four         | Pipeline Five         |
|----------------------------------------------|------------------|-----------------------|-----------------------|-----------------------|-----------------------|
| Raw Data QC                                  | Fastp            | Fastp                 | Fastp                 | Fastp                 | Fastp                 |
| Reads Alignment                              | STAR             | BWA                   | Bowtie2               | -                     | -                     |
| Reads counting                               | CIRCexplorer2    | CIRI                  | Find_circ             | Mapsplice             | Segemehl              |
| Data Processing (in house script)            | Python,Java,R    | Python,Java,R         | Python,Java,R         | Python,Java,R         | Python,Java,R         |
| Differential expression                      | edgeR            | edgeR                 | edgeR                 | edgeR                 | edgeR                 |
| Summary Report                               | MultiQC          | MultiQC               | MultiQC               | MultiQC               | MultiQC               |


### Dependencies
* Softwares
    * [Fastp](https://github.com/OpenGene/fastp)
    * [STAR](https://github.com/alexdobin/STAR)
    * [CIRCexplorer2](https://github.com/YangLab/CIRCexplorer2)
    * [BWA](https://github.com/lh3/bwa)
    * [CIRI](http://sourceforge.net/projects/ciri)
    * [Bowtie2](https://github.com/BenLangmead/bowtie2)
    * [Find_circ](https://github.com/marvin-jens/find_circ)
    * [Mapsplice](http://www.netlab.uky.edu/p/bioinfo/MapSplice2)
    * [Segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/)
    * [MultiQC](https://github.com/ewels/MultiQC)
    * Several R packages for downstream analysis.

