FROM ubuntu:24.04

LABEL version="0.0.8"
LABEL description="Conditional reciprocal best hits (CRBHits) in R"
LABEL maintainer="ullrich@evolbio.mpg.de"

RUN apt-get update && apt-get install -y adduser \
        && useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& adduser docker staff
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre \
                software-properties-common \
                dirmngr \
                gnupg \
                curl \
                ca-certificates \
	&& rm -rf /var/lib/apt/lists/*
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" | tee /etc/apt/sources.list.d/cran.list
RUN apt-get update \
    && apt-get install -y r-base=4.4.2-* \
    && apt-get clean
RUN R -e 'install.packages("BiocManager", repos="https://cloud.r-project.org/")' \
    && R -e 'BiocManager::install(version = "3.20")'
RUN apt-get install -y \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libglu1-mesa-dev \
        libgit2-dev \
        pandoc \
        libssh2-1-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev
RUN apt-get install -y \
        libpython3-dev \
        python3-dev \
        python3-pip \
        python3-venv
RUN R -e "install.packages('devtools')"
RUN R -e "message(BiocManager::version());BiocManager::install(version = \"3.20\");BiocManager::install(\"MSA2dist\");devtools::install_github(\"kullrich/CRBHits\", build_vignettes=FALSE, dependencies=TRUE); CRBHits::make_last(); CRBHits::make_dagchainer()"
RUN mkdir /usr/local/opt \
        && cd /usr/local/opt \
        && wget https://gitlab.com/mcfrith/last/-/archive/1595/last-1595.zip \
        && unzip last-1595.zip \
        && cd last-1595 \
        && make \
        && make install \
        && rm /usr/local/opt/last-1595.zip
RUN cd /usr/local/opt \
        && wget https://github.com/soedinglab/MMseqs2/releases/download/17-b804f/mmseqs-linux-avx2.tar.gz \
        && tar -xvf mmseqs-linux-avx2.tar.gz \
        && rm /usr/local/opt/mmseqs-linux-avx2.tar.gz
RUN cd /usr/local/opt \
        && wget https://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz \
        && tar -xvf diamond-linux64.tar.gz \
        && chmod +x diamond \
        && rm /usr/local/opt/diamond-linux64.tar.gz
CMD ["R"]
