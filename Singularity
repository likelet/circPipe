From:nfcore/base
Bootstrap:docker

%labels
    DESCRIPTION Singularity image containing all requirements for the nf-core/cirpipe pipeline
    VERSION 1.0dev

%environment
    PATH=/opt/conda/envs/nf-core-cirpipe-1.0dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
    conda install CIRCexplorer2
    conda install samtools
    conda install bowtie
    conda install bowtie2
    conda install fastp
    conda install star
    conda create -n tools_in_python3 python=3.6
    source activate tools_in_python3
    conda install mapsplice
    conda install segemehl
    source deactivate
    conda create -n tools_in_python2 python=2.7
    pip install pysam
    pip install numpy
    git clone http://github.com/marvin-jens/find_circ.git
    wget http://sourceforge.net/projects/ciri/files/CIRI-full/CIRI-full_v2.0.zip
    unzip CIRI-full_v2.0.zip
    export PATH=$PATH:/home/wqj/miniconda3/bin
    export PATH=$PATH:/home/wqj/miniconda3/envs/tools_in_python3/bin
    export PATH=$PATH:/home/wqj/miniconda3/envs/tools_in_python2/bin



