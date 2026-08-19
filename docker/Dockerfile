#FROM    public.ecr.aws/lts/ubuntu:latest
#FROM 	 rocker/tidyverse
#ENV     DEBIAN_FRONTEND noninteractive
FROM 	bioconductor/bioconductor_docker:RELEASE_3_22

## Otherwise timedatectl will get called which leads to 'no systemd' inside Docker
ENV TZ UTC

## Install dependencies for Python and R (following the official Rocker docs)
RUN apt-get update \
		&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre \
		python3-pip libxml2-dev libssl-dev libcurl4-openssl-dev \
		libv8-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
		libpng-dev libtiff5-dev libjpeg-dev software-properties-common dirmngr \
		&& apt-get clean \
		&& apt-get autoremove -y \
		&& apt-get autoclean -y 

# Install the required Python packages and add the CaTCH package
RUN     pip install pandas numpy biopython umi_tools
RUN     mkdir -p /packages/CaTCH && mkdir -p /tools/{data,scripts}
COPY    CaTCH /packages/CaTCH
COPY    data /tools/data
COPY    scripts /tools/scripts
ENV     PYTHONPATH /packages/
RUN     chmod +x /tools/scripts/python/*.py
ENV     PATH $PATH:/tools/scripts/python/

# Install the required R packages

ENV 	R_REMOTES_UPGRADE="never"

RUN		R -e "source('/tools/scripts/R/setup_R_env.R')"


# Install STAR mapper for STARsolo mapping option
ARG STAR_VERSION=2.7.11a

ENV PACKAGES build-essential g++ wget zlib1g-dev unzip default-jre

RUN set -ex

RUN apt-get install -y --no-install-recommends ${PACKAGES} && \
    	apt-get clean && \
		rm -rf /var/lib/apt/lists/* \
    	g++ --version && \
    	cd /home && \
    	wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.zip && \
    	unzip ${STAR_VERSION}.zip && \
    	cd STAR-${STAR_VERSION}/source && \
    	make STARstatic && \
    	cp STAR /usr/local/bin && \
    	cd /home && \
    	'rm' -rf STAR-${STAR_VERSION}

ENV PATH /home/.local/bin:${PATH}              

# Install Samtools
ARG SAMTOOLS_VERSION=1.19.1

ENV SAMTOOLS_INSTALL_DIR=/home/.local/bin/samtools

RUN cd /home && \
		wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
		tar --bzip2 -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
		cd samtools-${SAMTOOLS_VERSION} && \
		./configure --prefix=${SAMTOOLS_INSTALL_DIR} && \
		make && \
		make install && \
		cp ${SAMTOOLS_INSTALL_DIR}/bin/samtools /usr/local/bin && \
		cd /home && \
		rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2

# Install FASTQC
ARG FASTQC_VERSION=0.12.1
RUN cd /home && \
	wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip && \
	unzip -d /usr/local/lib/ fastqc_v${FASTQC_VERSION}.zip && \
	ln -s /usr/local/lib/FastQC/fastqc /usr/local/bin/. &&\
	cd /home && \
	rm -rf fastqc_v${FASTQC_VERSION}.zip

# Install MultiQC
RUN     pip install multiqc

LABEL   maintainer="sergej.nowoshilow@boehringer-ingelheim.com; joerg.fallmann@imp.ac.at"
LABEL   version="3.0"
LABEL   release-date="2023-10-04"