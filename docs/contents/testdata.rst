About Test data 
====================================

NOTE: A full packaged test data set can be downloaded from
http://cancerbio.info/pub/circpipe/testdata.tar.gz

Test data content
-----------------

::

    ├── compare.file
    ├── design.file
    ├── Fastq
    │   ├── Normal_S1_1.fq.gz
    │   ├── Normal_S1_2.fq.gz
    │   ├── Normal_S2_1.fq.gz
    │   ├── Normal_S2_2.fq.gz
    │   ├── Normal_S3_1.fq.gz
    │   ├── Normal_S3_2.fq.gz
    │   ├── Normal_S4_1.fq.gz
    │   ├── Normal_S4_2.fq.gz
    │   ├── Normal_S5_1.fq.gz
    │   ├── Normal_S5_2.fq.gz
    │   ├── Tumor_S1_1.fq.gz
    │   ├── Tumor_S1_2.fq.gz
    │   ├── Tumor_S2_1.fq.gz
    │   ├── Tumor_S2_2.fq.gz
    │   ├── Tumor_S3_1.fq.gz
    │   ├── Tumor_S3_2.fq.gz
    │   ├── Tumor_S4_1.fq.gz
    │   ├── Tumor_S4_2.fq.gz
    │   ├── Tumor_S5_1.fq.gz
    │   └── Tumor_S5_2.fq.gz
    ├── gencode.rsem.fpkm_m6Astatus_11_29.mat
    └── Genome
        ├── chr2.fa
        ├── gencode_chr2.v25.annotation.bed
        ├── gencode_chr2.v25.annotation.gtf
        └── hg38_chr2_gencode.txt

-  ``design.file`` store the exprimental design for performing
   comparision

::

                    Sample_id   Type
                    Normal_S1   Normal
                    Normal_S2   Normal
                    Normal_S3   Normal
                    Normal_S4   Normal
                    Normal_S5   Normal
                    Tumor_S1    Tumor
                    Tumor_S2    Tumor
                    Tumor_S3    Tumor
                    Tumor_S4    Tumor
                    Tumor_S5    Tumor

-  ``compare.file`` store the DE analysis comparision order

::

   Tumor_vs_Normal

-  ``gencode.rsem.fpkm_m6Astatus_11_29.mat`` is the mRNA expression
   matrix

-  ``Fastq`` folder contains paired, compressed fastq files, also known
   as raw reads.

-  ``Genome`` folder contains several files explained below:

   -  ``chr22.fa`` genome reference of chromosome 22 with fasta format.

   -  ``gencode_chr22.v25.annotation.gtf`` GTP files grepped from
      `GENCODE <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz>`__

   -  ``gencode_chr22.v25.annotation.bed`` GTP files grepped from
      `GENCODE <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz>`__

   -  ``hg38_chr2_gencode.txt`` annotation files of chromosome 22
