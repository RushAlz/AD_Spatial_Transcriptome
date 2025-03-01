FROM nvidia/cuda:10.2-runtime-ubuntu18.04
MAINTAINER Ricardo A. Vialle <ricardo_a_vialle@rush.edu>

# git clone https://github.com/snakemake/snakemake.git
ADD . /tmp/repo
WORKDIR /tmp/repo

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 curl ca-certificates gnupg2 squashfs-tools git g++

# Install Miniconda
RUN /bin/bash -c "curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh > mambaforge.sh && \
    bash mambaforge.sh -b -p /opt/conda && \
    rm mambaforge.sh"
    
RUN /bin/bash -c "mamba create -q -y -c conda-forge -c bioconda -n snakemake snakemake snakemake-minimal --only-deps && \
    source activate snakemake && \
    mamba install -q -y -c conda-forge singularity mkl-service && \
    conda clean --all -y && \
    which python && \
    pip install .[reports,messaging,google-cloud]"
RUN echo "source activate snakemake" > ~/.bashrc

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]

RUN pip3 install git+https://github.com/BayraktarLab/cell2location.git#egg=cell2location
RUN pip3 install anndata

ENV PATH /opt/conda/envs/snakemake/bin:${PATH}
