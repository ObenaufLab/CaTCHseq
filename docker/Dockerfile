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

# Install the required R packages
RUN     R -e "source('/tools/scripts/R/setup_R_env.R')"
              
LABEL   maintainer="sergej.nowoshilow@boehringer-ingelheim.com"
LABEL   version="2.2"
LABEL   release-date="2023-10-04"