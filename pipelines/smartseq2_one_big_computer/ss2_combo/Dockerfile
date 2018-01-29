FROM ubuntu:17.10

ENV version 2.10.10
WORKDIR /usr/picard
ADD https://github.com/broadinstitute/picard/releases/download/${version}/picard.jar ./picard.jar

# Install R dependencies for Picard
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y r-base

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('reshape')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('ggplot2')"

RUN mkdir -p /opt/tools/
WORKDIR /opt/tools

RUN \
 apt update && \ 
 apt install -y \
  liblzma-dev \
  libbz2-dev \
  cmake \
  automake \
  curl \
  libboost-all-dev \
  wget \
  build-essential \    
  gcc-multilib \
  zlib1g-dev \
  libxml2-dev \
  libncurses5-dev \
  zip unzip \
  samtools \
  git \
  r-base \
  r-base-core \
  r-base-dev
  
RUN \
  wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip && \
  unzip hisat2-2.1.0-source.zip && \
  cd hisat2-2.1.0 && \
  make && \
  cp hisat2* /usr/local/bin

COPY  ./make_grch38_snp_tran_*.sh hisat2-2.1.0/
# set ENV PATH, run hisat2 by type 'HISAT2' from cmd
ENV PATH="/opt/tools/hisat2-2.1.0:${PATH}"

# gffread is gtf/gtt tools 
WORKDIR /opt/tools/gffread
RUN \
  git clone https://github.com/gpertea/gclib && \
  git clone https://github.com/gpertea/gffread && \
  cd gffread && \
  make && \ 
  cp gffread /usr/local/bin/

RUN apt-get update --fix-missing && \
  apt-get install -y \
  build-essential \
  gcc-multilib \
  apt-utils \
  zlib1g-dev \
  libxml2-dev \
  curl \
  wget \
  git \
  perl \
  perl-base \
  libbz2-dev \
  cmake automake \
  libboost-all-dev \
  libncurses5-dev \
  r-base \
  r-base-core \
  r-base-dev \
  libncurses5-dev \
  bowtie \
  bowtie2 

# Compile and install STAR
WORKDIR /usr/local
RUN wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz && \
    tar -xf 2.5.3a.tar.gz
WORKDIR STAR-2.5.3a/bin/Linux_x86_64_static
RUN cp STAR /usr/local/bin

# Install RSEM 
WORKDIR /usr/local/
RUN git clone https://github.com/deweylab/RSEM.git
WORKDIR /usr/local/RSEM
RUN git checkout v1.3.0
RUN make 
RUN make ebseq
ENV PATH /usr/local/RSEM:$PATH

RUN apt-get install -y openjdk-8-jre-headless
WORKDIR /opt/tools
