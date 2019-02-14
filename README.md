# circPipe
**circRNA analysis pipe**

[![Build Status](https://travis-ci.org/likelet/cirpipe.svg?branch=master)](https://travis-ci.org/likelet/cirpipe)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![Docker](https://img.shields.io/docker/automated/likelet/cirpipe.svg)](https://hub.docker.com/r/likelet/cirpipe)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)
[![Documentation Status](https://readthedocs.org/projects/circpipe/badge/?version=latest)](https://circpipe.readthedocs.io/en/latest/?badge=latest)


### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation   

A full tutorial of CircPipe can be found at Wiki page of this project. plz go to the https://github.com/likelet/circPipe/wiki


### Pipeline Steps

CircPipe allows you to run arbitrary pipelines among five pipelines.
Choose between workflows by using `--selectTools` or not(default) .

| Step                                         | Pipeline One     | Pipeline Two          | Pipeline Three        | Pipeline Four         | Pipeline Five         |Pipeline Six         |
|----------------------------------------------|------------------|-----------------------|-----------------------|-----------------------|-----------------------|---------------------|
| Raw Data QC                                  | Fastp            | Fastp                 | Fastp                 | Fastp                 | Fastp                 |Fastp                |
| Reads Alignment                              | STAR             | BWA                   | Bowtie2               | -                     | -                     | -                   |
| Reads counting                               | CIRCexplorer2    | CIRI                  | Find_circ             | Mapsplice             | Segemehl              |KNIFE                |
| Data Processing (in house script)            | Python,Java,R    | Python,Java,R         | Python,Java,R         | Python,Java,R         | Python,Java,R         |Python,Java,R        |
| Differential expression                      | edgeR            | edgeR                 | edgeR                 | edgeR                 | edgeR                 |edgeR                |
| Summary Report                               | MultiQC          | MultiQC               | MultiQC               | MultiQC               | MultiQC               |MultiQC              |


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
    * [KNIFE](https://github.com/lindaszabo/KNIFE)
    * [MultiQC](https://github.com/ewels/MultiQC)
    * Several R packages for downstream analysis.

