# nf-core/cirpipe
**cirRNA analysis pipe**

[![Build Status](https://travis-ci.org/nf-core/cirpipe.svg?branch=master)](https://travis-ci.org/nf-core/cirpipe)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/cirpipe.svg)](https://hub.docker.com/r/nfcore/cirpipe)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The nf-core/cirpipe pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)


### Pipeline Steps

The pipeline allows you to choose between running either replicates or without replicates.
Choose between workflows by using `--without_replicate` or not(default) .

| Step                                         | Pipeline One     | Pipeline Two          | Pipeline Three        | Pipeline Four         | Pipeline Five         | Pipeline Six          |
|----------------------------------------------|------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|
| Raw Data QC                                  | Fastp            | Fastp                 | Fastp                 | Fastp                 | Fastp                 | Fastp                 |
| Reads Alignment                              | STAR             | BWA                   | Bowtie2               | -                     | -                     | Bowtie2               |
| Reads counting                               | CIRCexplorer2    | CIRI                  | Find_circ             | Mapsplice             | Segemehl              | AutoCirc              |
| Data Processing (in house script)            | Python,JAVA,R    | Python,JAVA,R         | Python,JAVA,R         | Python,JAVA,R         | Python,JAVA,R         | Python,JAVA,R         |
| Differential expression                      | edgeR            | edgeR                 | edgeR                 | edgeR                 | edgeR                 | edgeR                 |
| Summary Report                               | MultiQC          | MultiQC               | MultiQC               | MultiQC               | MultiQC               | MultiQC               |


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
    * [AutoCirc](https://github.com/chanzhou/AutoCirc)
    * [MultiQC](https://github.com/ewels/MultiQC)
    * Several R packages for downstream analysis.