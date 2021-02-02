FROM conda/miniconda3

LABEL maintainer="andreas.hoek@bio.uni-giessen.de"

RUN conda install -y -c conda-forge mamba
RUN mamba create -y -c conda-forge -c bioconda -c anaconda -n WASP snakemake=5.20.1 fastqc star subread samtools umi_tools simplejson python-levenshtein py more-itertools zip

COPY Snakefile /WASP/Snakefile
COPY wasp_docker.sh /WASP/wasp_docker.sh
COPY Scripts /WASP/Scripts/

WORKDIR /WASP/
ENTRYPOINT ["/bin/bash", "wasp_docker.sh", "-p"]
