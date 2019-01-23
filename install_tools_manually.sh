#!/usr/bin/sh

# install CIRI
echo "Downloading CIRI..."
wget http://sourceforge.net/projects/ciri/files/CIRI-full/CIRI-full_v2.0.zip 
echo "Done!"
echo "Installing CIRI..."
unzip CIRI-full_v2.0.zip
rm CIRI-full_v2.0.zip
mv CIRI-full_v2.0 CIRI
sed -i "1i\\#\!/usr/bin/perl" ./CIRI/bin/CIRI_v2.0.6/CIRI2.pl
chmod a+x ./CIRI/bin/CIRI_v2.0.6/*
echo "export PATH=\"$(pwd)/CIRI/bin/CIRI_v2.0.6:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
echo "Finish CIRI installation!"

# install Find_circ
echo "Downloading Find_circ..."
git clone http://github.com/marvin-jens/find_circ.git 
echo "Done!"
echo "Installing Find_circ..."
sed -i '1d' ./find_circ/unmapped2anchors.py
sed -i "1i\\#\!/usr/bin/python" ./find_circ/unmapped2anchors.py
sed -i '1d' ./find_circ/find_circ.py
sed -i "1i\\#\!/usr/bin/python" ./find_circ/find_circ.py
chmod a+x ./find_circ/*
echo "export PATH=\"$(pwd)/find_circ:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
echo "Finish Find_circ installation!"

#install the tools by miniconda....
echo "Make sure you have install miniconda before!..."
echo "Installing the tools by miniconda...."
conda install fastp=0.19.5
conda install hisat2=2.1.0
conda install htseq=0.11.2
conda install multiqc=1.6
conda install pysam=0.15.1
conda install perl=5.26.2.1
conda install numpy=1.15.2
conda install segemehl=0.2.0
conda install mapsplice=2.2.0
conda install bwa=0.7.17
conda install nextflow=0.32.0
conda install graphviz=2.38.0
conda install bowtie=1.2.2
conda install bowtie2=2.3.4.3
conda install samtools=1.9
conda install star=2.6.1b
conda install circexplorer2=2.3.3
conda install r-ggplot2=3.1.0
conda install r-pheatmap=1.0.10
conda install bioconductor-edger=3.22.5
conda install bioconductor-deseq2=1.20.0
conda install bioconductor-limma=3.36.5
conda install bioconductor-chipseeker=1.16.1
conda install bioconductor-fgsea=1.6.0
conda install bioconductor-bsgenome.hsapiens.ucsc.hg38=1.4.1
conda install bioconductor-clusterprofiler=3.8.1
conda install bioconductor-reactomepa=1.24.0
conda install bioconductor-pathview=1.20.0
conda install r-dplyr=0.7.6
conda install r-readr=1.1.1
conda install r-rmarkdown=1.10
conda install r-ggsci=2.9
conda install r-circlize=0.4.5
conda install r-rcolorbrewer=1.1_2
conda install r-venndiagram=1.6.20
conda install bedtools=2.27.1
conda install biopython=1.72
conda install bcbiogff=0.6.4
conda install perl=5.26.2.1
conda install r-markdown=0.9
conda install scipy=1.2.0

echo "Finish all the installation!"
