FROM    public.ecr.aws/lts/ubuntu:latest

ENV     DEBIAN_FRONTEND noninteractive

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
	&& rm -rf /var/lib/apt/lists/*

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Install R in pinned version
ENV R_BASE_VERSION 4.3.2

## Now install R and littler, and create a link for littler in /usr/local/bin
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
			libopenblas0-pthread \
	littler \
			r-cran-docopt \
			r-cran-littler \
	r-base=${R_BASE_VERSION}-* \
	r-base-dev=${R_BASE_VERSION}-* \
			r-base-core=${R_BASE_VERSION}-* \
	r-recommended=${R_BASE_VERSION}-* \
	&& chown root:staff "/usr/local/lib/R/site-library" \
	&& chmod g+ws "/usr/local/lib/R/site-library" \
	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installDeps.r /usr/local/bin/installDeps.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

RUN     wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

RUN     add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN     apt install r-base r-base-dev -y

# Install the required Python packages and add the CaTCH package
RUN     pip install pandas numpy biopython
RUN     mkdir -p /packages/CaTCH && mkdir -p /tools/{data,scripts}
COPY    CaTCH /packages/CaTCH
COPY    data /tools/data
COPY    scripts /tools/scripts
ENV     PYTHONPATH /packages/
RUN     chmod +x /tools/scripts/python/*.py
ENV     PATH $PATH:/tools/scripts/python/

# Install STAR mapper for STARsolo mapping option
ARG STAR_VERSION=2.7.11a

ENV PACKAGES gcc g++ make wget zlib1g-dev unzip vim

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

# Install the required R packages
RUN     R -e "source('/tools/scripts/R/setup_R_env.R')"
              
LABEL   maintainer="sergej.nowoshilow@boehringer-ingelheim.com"
LABEL   version="3.0"
LABEL   release-date="2023-10-04"