Home 
=============================

In a nutshell, `circPipe <https://github.com/likelet/circPipe>`_
(circRNAs Pipeline) aims at data exploration of circRNAs. It begins with
the raw sequencing data and then following a step of quality control. We
recommend our users to use the paired-end sequences and the sequencing
depth should be more than 70M.The length of the reads should be longer
than 50bp, and longer than 100bp is the best.We absorb five kinds of
common software/work to detect circRNAs, including
`Circexplorer2 <https://circexplorer2.readthedocs.io/en/latest/>`_,
`CIRI <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0571-3>`_,
`Find_circ <https://github.com/marvin-jens/find_circ>`_,
`Mapsplice <http://www.netlab.uky.edu/p/bioinfo/MapSplice2>`_,
and `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_. Users can
choose one, several or all appropriate software according to their
project. If our users chose several kinds of software, our pipeline
will then combine all these results in to a sole matrix. To interpret the data, we design three modules to
explore the identified circRNAs as well as protein coding ones,
including conventional annotation, differential circRNA expression
analysis, annotation and pathway enrichment. Plots and tables of
analysis module are presented in a HTML file via
`Rmarkdown <https://rmarkdown.rstudio.com/>`__.

More information can be found in the The cirPipe workflow page, or in
the project GitHub wiki. Please be sure that you have all dependencies
software or tools preinstalled in you system. Otherwise, we recommended
that users employ docker or singularity containers to run the pipe.

This wiki includes several tutorials, plz following the step by step
tutorial to run cirPipe in your local server or cluster. This tutorial
explain how to set the parameters in the ``nextflow.config`` file, and
describe the files that will be produced in output, while at this page
you can know more about How to read the logs.


Workflows
~~~~~~~~~

Our Pipeline Steps are showed in the below chart. CircPipe allows you to
run arbitrary pipelines among five pipelines. Choose between workflows
by using ``--selectTools``.

+--------------------+----------+-----------+-----------+-----------+-----------+
| Step               | Pipelin  | Pipeline  | Pipeline  | Pipeline  | Pipeline  |
|                    | e        | Two       | Three     | Four      | Five      |
|                    | One      |           |           |           |           |
+====================+==========+===========+===========+===========+===========+
| Raw Data QC        | Fastp    | Fastp     | Fastp     | Fastp     | Fastp     |
+--------------------+----------+-----------+-----------+-----------+-----------+
| Reads Alignment    | STAR     | BWA       | Bowtie2   | -         | -         |
+--------------------+----------+-----------+-----------+-----------+-----------+
| cirRNA             | CIRCexp  | CIRI      | Find\_cir | Mapsplice | Segemehl  |
| identification     | lorer2   |           | c         |           |           |
+--------------------+----------+-----------+-----------+-----------+-----------+
| Data Processing    | Python,  | Python,Ja | Python,Ja | Python,Ja | Python,Ja |
| (in house script)  | Java,R   | va,R      | va,R      | va,R      | va,R      |
+--------------------+----------+-----------+-----------+-----------+-----------+
| Differential       | edgeR    | edgeR     | edgeR     | edgeR     | edgeR     |
| expression         |          |           |           |           |           |
+--------------------+----------+-----------+-----------+-----------+-----------+
| Summary Report     | MultiQC  | MultiQC+R | MultiQC+R | MultiQC+R | MultiQC+R |
|                    | +R       | markdown  | markdown  | markdown  | markdown  |
|                    | markdow  |           |           |           |           |
|                    | n        |           |           |           |           |
+--------------------+----------+-----------+-----------+-----------+-----------+
