FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/circpipe pipeline"

COPY environment.yml environment1.yml environment2.yml ./

ENV PATH /opt/conda/bin:$PATH
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/tools_in_python2/bin:$PATH
RUN conda env create -f /environment1.yml -n tools_in_python2 python=2.7 && conda clean -a

ENV PATH /opt/conda/envs/tools_in_python3/bin:$PATH
RUN conda env create -f /environment2.yml -n tools_in_python3 python=3.6 && conda clean -a

ENV PATH /opt/conda/envs/nf-core-cirpipe-1.0dev/bin:$PATH

#install ciri
RUN wget http://sourceforge.net/projects/ciri/files/CIRI-full/CIRI-full_v2.0.zip && \
    unzip CIRI-full_v2.0.zip && \
    rm CIRI-full_v2.0.zip && \
    mv CIRI-full_v2.0 CIRI

#install find_circ
RUN git clone http://github.com/marvin-jens/find_circ.git

#install KNIFE
RUN wget https://github.com/lindaszabo/KNIFE/archive/v1.4.tar.gz && \
    tar zxvf v1.4.tar.gz && \
    rm v1.4.tar.gz && \
    mv KNIFE-1.4 KNIFE && \
    cd KNIFE/circularRNApipeline_Standalone/analysis && \
    chmod a+x *

