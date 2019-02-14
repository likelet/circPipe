Prepare Input files
===================

To use this pipeline, users should prepare the raw read sequences with
or without biological replicate in fastq format first. In addition, we
recommend our users to use the paired-end sequences and the sequencing
depth should be more than 70M. The length of the reads should be longer
than 50bp, and longer than 100bp is the best.

To use the five kinds of software and annotation, users should choose an
appropriate reference index as one of the input files. If users don't
provide a specific index file, the pipeline will use Hg19 index file by
default.

For the step of differential gene expression analysis, users should
provide ``design.txt`` as a design file and ``compare.txt`` as a compare
file.

The parameters of each input file are described in the section of
``Parameters`` in detail, users can refer to that section and choose the
input files as their specific needs.

Input file
~~~~~~~~~~

-  | ``design.txt``
   | sampleInfor presents the experimental design of your data set, it
     is just like a design file of ``DESeq2`` and ``EdgeR`` input.

   ::

       Sample  Type
       P1003NA N
       P1003TA T
       P1162NA N
       P1162TA T
       P1408NA N
       P1408TA T
       P1527NA N

-  ``compare.txt`` specify which group to compare in your differential
   expression analysis

   ::

       T_vs_N

   ``T`` and ``N`` are the identical strings as the ``Type`` column in
   ``design.txt``.
