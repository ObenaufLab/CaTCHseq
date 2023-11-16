#FROM    public.ecr.aws/lts/ubuntu:latest
#FROM 	 rocker/tidyverse
#ENV     DEBIAN_FRONTEND noninteractive
FROM 	bioconductor/bioconductor_docker:RELEASE_3_18

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
		&& apt-get autoclean -y \
		&& rm -rf /var/lib/apt/lists/*

# Install the required Python packages and add the CaTCH package
RUN     pip install pandas numpy biopython
RUN     mkdir -p /packages/CaTCH && mkdir -p /tools/{data,scripts}
COPY    CaTCH /packages/CaTCH
COPY    data /tools/data
COPY    scripts /tools/scripts
ENV     PYTHONPATH /packages/
RUN     chmod +x /tools/scripts/python/*.py
ENV     PATH $PATH:/tools/scripts/python/

# Install the required R packages
			
RUN		R -e "source('/tools/scripts/R/setup_R_env.R')"


# Install STAR mapper for STARsolo mapping option
ARG STAR_VERSION=2.7.11a

ENV PACKAGES build-essential g++ wget zlib1g-dev unzip

RUN set -ex

RUN apt-get install -y --no-install-recommends ${PACKAGES} && \
    	apt-get clean && \
    	g++ --version && \
    	cd /home && \
    	wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.zip && \
    	unzip ${STAR_VERSION}.zip && \
    	cd STAR-${STAR_VERSION}/source && \
    	make STARstatic && \
    	mkdir -p /home/.local/bin && \
    	cp STAR /home/.local/bin && \
    	cd /home && \
    	'rm' -rf STAR-${STAR_VERSION}

ENV PATH /home/.local/bin:${PATH}              

LABEL   maintainer="sergej.nowoshilow@boehringer-ingelheim.com; joerg.fallmann@imp.ac.at"
LABEL   version="3.0"
LABEL   release-date="2023-10-04"