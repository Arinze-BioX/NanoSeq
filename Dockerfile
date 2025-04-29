FROM  ubuntu:18.04 as builder

USER  root

# ALL tool versions used by opt-build.sh
ENV VER_SAMTOOLS="1.18"
ENV VER_HTSLIB="1.18"
ENV VER_BCFTOOLS="1.18"
ENV VER_VERIFYBAMID="2.0.1"
ENV VER_LIBDEFLATE="v1.18"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends locales
RUN apt-get install -yq --no-install-recommends g++
RUN apt-get install -yq --no-install-recommends ca-certificates
RUN apt-get install -yq --no-install-recommends wget

# install latest cmake so opt-build.sh works - the initial installs will also help install R
RUN apt-get install -yq --no-install-recommends software-properties-common lsb-release
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
RUN apt-add-repository "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main"
RUN apt-get install -yq --no-install-recommends cmake=3.25.2-0kitware1ubuntu18.04.1

RUN apt-get install -yq --no-install-recommends make
RUN apt-get install -yq --no-install-recommends pkg-config

# if ubuntu 18.04
RUN apt install -yq --no-install-recommends dirmngr
RUN apt-get install -yq --no-install-recommends zlib1g-dev
RUN apt-get install -yq --no-install-recommends libbz2-dev
RUN apt-get install -yq --no-install-recommends liblzma-dev
RUN apt-get install -yq --no-install-recommends libcurl4-openssl-dev
RUN apt-get install -yq --no-install-recommends libncurses5-dev
RUN apt-get install -yq --no-install-recommends libssl-dev
RUN apt-get install -yq --no-install-recommends libblas-dev
RUN apt-get install -yq --no-install-recommends liblapack-dev
RUN apt-get install -yq --no-install-recommends gfortran
RUN apt-get install -yq --no-install-recommends libxml2-dev
RUN apt-get install -yq --no-install-recommends libgsl-dev
RUN apt-get install -yq --no-install-recommends libperl-dev
RUN apt-get install -yq --no-install-recommends libpng-dev


RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# build tools from other repos
ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

# build the tools in this repo, separate to reduce build time on errors
COPY . .
RUN bash build/opt-build-local.sh $OPT

FROM ubuntu:18.04

LABEL maintainer="okafor.ae@gmail.com" \
      version="1.0.1" \
      description="nanorateseq docker"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
curl \
wget \
make \
g++ \
gcc \
gfortran \
libblas-dev \
liblapack-dev \
ca-certificates \
time \
zlib1g \
libz-dev \
libxml2 \
libgsl23 \
libperl5.26 \
libcapture-tiny-perl \
libfile-which-perl \
libpng16-16 \
parallel \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN apt install -yq --no-install-recommends software-properties-common dirmngr

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

# install miniforge
ENV HOME /home/ubuntu
ENV PATH /home/ubuntu/miniforge3/bin:${PATH}
ENV CONDA_DIR $HOME/miniforge3
RUN wget https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-$(uname -m).sh
RUN bash Miniforge3-24.11.3-2-Linux-$(uname -m).sh -b -p $(pwd)/miniforge3
RUN rm Miniforge3-24.11.3-2-Linux-$(uname -m).sh

# initialize conda and install mamba
RUN $HOME/miniforge3/bin/conda init
SHELL ["/bin/bash", "-i", "-c"] 
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -y mamba=2.0.7

# Install other packages from conda using mamba
RUN mamba install -y -c conda-forge -y r-base=4.3.3 r-epitools=0.5_10.1 r-ggplot2=3.5.1 r-data.table=1.17.0 r-gridextra=2.3 perl=5.32.1 perl-file-which=1.24 r-upsetr=1.4.0 r-ggplot2=3.5.1
RUN mamba install -c bioconda -y bioconductor-deepsnv=1.48.0 r-vcfr=1.15.0 r-seqinr=4.2_36 snakemake=9.1.1 snakemake-executor-plugin-slurm bwa=0.7.19 biobambam=2.0.185 verifybamid2=2.0.1 gatk=3.8
RUN mamba install -y r::r-vgam=1.1_9
RUN mamba install -y -c bioconda perl-getopt-long=2.58 perl-pod-usage=2.05 perl-capture-tiny=0.48

ENV PATH /home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${PATH}

CMD ["/bin/bash", "-i"]