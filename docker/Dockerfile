FROM    public.ecr.aws/lts/ubuntu:latest

RUN     apt-get update && \
        apt-get upgrade -y

ENV     DEBIAN_FRONTEND noninteractive
RUN     apt-get install -y python3-pip r-base libxml2-dev libssl-dev libcurl4-openssl-dev \
        libv8-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
        libpng-dev libtiff5-dev libjpeg-dev

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