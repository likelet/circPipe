Parameters
----------

.. warning :: Those parameters would cover the setting from ``nextflow.config`` file .

Mandatory
^^^^^^^^^

.. note:: plz configure those options in ``nextflow.config`` or ``docker.config`` file .

+--------------+------------------+----------------+
| Name         | Example/Default  | Description    |
|              | value            |                |
+==============+==================+================+
| --reads      | ``./*{1,2}.fq.gz | input raw      |
|              | ``               | paired reads   |
+--------------+------------------+----------------+
| --designfile | ``FALSE``        | A txt file     |
|              |                  | that stored    |
|              |                  | experimental   |
|              |                  | design         |
|              |                  | information,   |
|              |                  | plz see        |
|              |                  | details from   |
|              |                  | ``--designfile |
|              |                  | ``             |
|              |                  | section above  |
+--------------+------------------+----------------+
| --comparefile| ``FALSE``        | A txt file     |
|              |                  | that stored    |
|              |                  | experimental   |
|              |                  | compare        |
|              |                  | information,   |
|              |                  | plz see        |
|              |                  | details from   |
|              |                  | ``--comparefil |
|              |                  | e``            |
|              |                  | section above  |
+--------------+------------------+----------------+

Configuration
^^^^^^^^^^^^^

    (paths to references, softwares and special environments. Only need
    to be set at the first time you run the pipeline) .

+--------------+-----------------+----------------+
| Name         | Example/Default | Description    |
|              | value           |                |
+==============+=================+================+
| --genomefile | ``path/to/refdi | Path to Fasta  |
|              | r/chr2.fa``     | reference(requ |
|              |                 | ired           |
|              |                 | if not set in  |
|              |                 | config file)   |
+--------------+-----------------+----------------+
| --gtffile/-- | ``path/to/refdi | Different      |
| bedfile/--an | r/hg38_gencode. | annotation     |
| notationfile | txt``           | files from     |
|              |                 | GENCODE        |
|              |                 | database for   |
|              |                 | annotating     |
|              |                 | circRNAs. e.g. |
|              |                 | [gencode.v25.a |
|              |                 | nnotation.gtf] |
|              |                 | /[gencode.v25. |
|              |                 | annotation.bed |
|              |                 | ]/[hg38\_genco |
|              |                 | de.txt]        |
+--------------+-----------------+----------------+
| --ciridir/-- | ``/path/to/tool | Home folder of |
| find\_circdi | s/directory``   | ciri/find\_cir |
| r/--mapsdir/ |                 | c/mapsplice    |
|              |                 | installed      |
|              |                 | location       |
|              |                 |                |
+--------------+-----------------+----------------+
|--genomebuild |``'hg19'``       |specific genome |
|              | 'GRCh38', 'hg10'| build for      |
|              | available       |circplot        |
+--------------+-----------------+----------------+

Optional
^^^^^^^^

+--------------+-----------------+----------------+
| Name         | Default value   | Description    |
+==============+=================+================+
|-profile      | ``standard``    | Configuration  |
|              |                 | profile to use.|
|              |                 | Available:     |
|              |                 | standard,      |
|              |                 | conda, docker, |
|              |                 | singularity,   |
|              |                 | test           |
+--------------+-----------------+----------------+
|--singleEnd   | ``false``       | specify that   |
|              |                 | the reads are  |
|              |                 | single ended   |
+--------------+-----------------+----------------+
|--selectTools | ``1``           | specify which  |
|              |                 | tools should   |
|              |                 | be use. ``1``  |
|              |                 | for            |
|              |                 | circexplorer2, |
|              |                 | ``2`` for      |
|              |                 | ciri, ``3``    |
|              |                 | for            |
|              |                 | find\_circ,    |
|              |                 | ``4`` for      |
|              |                 | mapsplice,     |
|              |                 | ``5`` for      |
|              |                 | segemehl,      |
|              |                 | For example,   |
|              |                 | you can set    |
|              |                 | ``1,2,3,4,5``  |
|              |                 | for running    |
|              |                 | five tools in  |
|              |                 | the same time. |
+--------------+-----------------+----------------+
| --skipDE     | ``false``       | skip           |
|              |                 | differential   |
|              |                 | expression     |
|              |                 | analysis       |
+--------------+-----------------+----------------+         
| --outdir     | ``./Result``    | the output     |
|              |                 | directory of   |
|              |                 | the results    |
+--------------+-----------------+----------------+
| --mRNA       | ``path/to/genco | Path to the    |
|              | de.rsem.fpkm_m6 | mRNA           |
|              | Astatus_11_29.m | expression     |
|              | at``            | matrix. Only   |
|              |                 | need to be set |
|              |                 | when you want  |
|              |                 | to do the      |
|              |                 | correlation.   |
+--------------+-----------------+----------------+

Detailed instruction of parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``--reads``

   suffix of your raw reads file. For example, ``*_{1,2}.fq.gz`` for
   paired end reads file ``sampleA_1.fq.gz`` and ``sampleA_2.fq.gz``

-  ``--outdir``

   path to the result output directory

-  ``--selectTools``

   specify which tools should be use. ``1`` for circexplorer2, ``2`` for
   ciri, ``3`` for find\_circ, ``4`` for mapsplice, ``5`` for segemehl.
   For example, you can set
   ``--selectTools='1,2,3,4,5'`` for running five tools in the same
   time.

-  ``--designfile``

   design file

-  ``--comparefile``

   compare file

-  ``--mRNA``

   mRNA expression matrix file

-  ``--gtffile``

   gtf file for building your STAR index, running CIRI and Mapsplice,
   running annotation. For example, ``gencode.v25.annotation.gtf``.

-  ``--genomefile``

   whole genome reference sequence in ``.fa`` format for running
   CIRCexplorer2, CIRI, Segemehl, Find\_circ. For example,
   ``genome.fa``.

-  ``--annotationfile``

   annotation file of genome in ``.txt`` format for running
   CIRCexplorer2. For example, ``hg38_gencode.txt``.

-  ``--singleEnd``

   ``true`` when using a single End reads input, default ``false``

Configure profiles 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As a nextflow-based analysis pipeline, CircPipe allow users to edit configure file ``nextflow.config`` to set the index files and default file path parameters instead of typing them into the command line.

To configure, please go to ``params`` line, and set the following information of various file locations and system environment settings.
Here, we use the test.config as an example.

.. code:: groovy

    params {

  //reads files
  reads = "Fastq/*_{1,2}.fq.gz"

  // design files and compare file 
  designfile="design.file"
  comparefile="compare.txt"

  currentPath="/data2/zhaoqi/circlePipetest/"
  //the necessary reference
  refmapsplice = false
  annotationfile = "${currentPath}Genome/hg19_chr2_refseq.txt" // for using circexplorer2, can be obtained from 
  genomefile = "${currentPath}Genome/hg19_chr2.fa"
  faifile = "${currentPath}Genome/hg19_chr2.fa.fai"
  gtffile = "${currentPath}Genome/hg19_chr2.gencode.annotation.gtf"
  mRNA = ""

  // index files for each software 
  starindex = "${currentPath}Genome/hg19_chr2_starindex"// path and prefix 
  segindex = "${currentPath}Genome/hg19_chr2" // path only 
  hisat2_index = "${currentPath}Genome/hg19_chr2_hisat2_index/hg19_chr2"
  bowtie2index = "${currentPath}Genome/hg19_chr2" // path and prefix 
  bowtieindex = "${currentPath}Genome/hg19_chr2" // path and prefix
  bwaindex = "${currentPath}Genome/hg19_chr2" //path and prefix
  skipDE = false
  
  }
