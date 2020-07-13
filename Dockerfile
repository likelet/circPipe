FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/circpipe pipeline"

COPY environment.yml ./

RUN conda create -n find_circ python=2.7

ENV PATH /opt/conda/bin:$PATH

RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/cirpipe-1.0dev/bin:$PATH

#install ciri
RUN wget http://sourceforge.net/projects/ciri/files/CIRI-full/CIRI-full_v2.0.zip && \
    unzip CIRI-full_v2.0.zip && \
    rm CIRI-full_v2.0.zip && \
    mv CIRI-full_v2.0 CIRI && \
    sed -i "1i\\#\!/usr/bin/perl" ./CIRI/bin/CIRI_v2.0.6/CIRI2.pl && \
    chmod a+x ./CIRI/bin/CIRI_v2.0.6/* && \
    echo "export PATH=\"$(pwd)/CIRI/bin/CIRI_v2.0.6:\$PATH\"" >> ~/.bashrc && \
    . ~/.bashrc

ENV PATH /CIRI/bin/CIRI_v2.0.6:$PATH

SHELL ["/bin/bash", "-c"]
RUN source activate find_circ && \
    pip install pysam && \
    conda install -c bioconda bowtie2 samtools

#install find_circ
RUN git clone http://github.com/marvin-jens/find_circ.git && \
    sed -i '1d' ./find_circ/unmapped2anchors.py && \
    sed -i "1i\\#\!/usr/bin/env python" ./find_circ/unmapped2anchors.py && \
    sed -i '1d' ./find_circ/find_circ.py && \
    sed -i "1i\\#\!/usr/bin/env python" ./find_circ/find_circ.py && \
    chmod a+x ./find_circ/* && \
    echo "export PATH=\"$(pwd)/find_circ:\$PATH\"" >> ~/.bashrc && \
    . ~/.bashrc

ENV PATH /find_circ:$PATH