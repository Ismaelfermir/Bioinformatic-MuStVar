FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y apt-transport-https
RUN apt-get install -y wget
RUN apt-get install -y build-essential
RUN apt-get install -y python2.7
RUN apt-get install -y sudo
RUN apt-get install -y unzip
RUN apt-get install -y tar
RUN apt-get install -y python-pip
RUN apt-get install -y libx11-dev
RUN apt-get install -y libxpm-dev
RUN apt-get install -y libxft-dev
RUN apt-get install -y libxext-dev
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y libncursesw5-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev
RUN apt-get install -y git binutils
RUN apt-get install -y g++
RUN apt-get install -y automake
RUN apt-get install -y git
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libopenblas-base
RUN apt-get install -y libopenblas-dev
RUN apt-get install -y cpanminus
RUN apt-get install -y fort77
RUN apt-get install -y xorg-dev
RUN apt-get install -y liblzma-dev  libblas-dev gfortran
RUN apt-get install -y gcc-multilib
RUN apt-get install -y gobjc++
RUN apt-get install -y aptitude
RUN aptitude install -y libreadline-dev
RUN apt-get install -y libcurl4-gnutls-dev
RUN apt-get install -y nano
RUN apt-get install -y curl
RUN apt-get install -y gcc
RUN apt-get install -y bzip2
RUN sudo apt-get install -y build-essential git libncurses-dev

ADD . /home

WORKDIR /home

RUN cd /home/tools && wget https://ftp.pcre.org/pub/pcre/pcre-8.41.zip
RUN cd /home/tools && unzip pcre-8.41.zip
RUN rm -r /home/tools/pcre-8.41.zip
RUN cd /home/tools/pcre-8.41 && ./configure --prefix=/usr --docdir=/usr/share/doc/pcre-8.42 --enable-unicode-properties --enable-pcre16 --enable-pcre32 --enable-pcregrep-libz --enable-pcregrep-libbz2 --enable-pcretest-libreadline --disable-static && make
RUN cd /home/tools/pcre-8.41 && make install && mv -v /usr/lib/libpcre.so.* /lib && ln -sfv ../../lib/$(readlink /usr/lib/libpcre.so) /usr/lib/libpcre.so

RUN cd /home/tools && wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
RUN tar -xvjf /home/tools/htslib-1.10.2.tar.bz2
RUN cd /home/tools/htslib-1.10.2 && ./configure
RUN cd /home/tools/htslib-1.10.2 && make
RUN cd /home/tools/htslib-1.10.2 && make install
RUN rm -r /home/tools/htslib-1.10.2.tar.bz2

RUN cd /home/tools && wget hhttps://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN tar -xvjf /home/tools/samtools-1.10.tar.bz2
RUN cd /home/tools/samtools-1.10 && ./configure
RUN cd /home/tools/samtools-1.10 && make
RUN rm -r /home/tools/samtools-1.10.tar.bz2
RUN cd /home/tools/samtools-1.10 && cp samtools /usr/bin

RUN cd /home/tools && wget https://github.com/broadinstitute/gatk/releases/download/4.0.6.0/gatk-4.0.6.0.zip
RUN cd /home/tools && unzip gatk-4.1.9.0.zip
RUN cd /home/tools && rm -r gatk-4.1.4.1.zip
RUN cd /home/tools/gatk-4.1.9.0 && cp gatk /usr/bin
RUN export GATK_LOCAL_JAR=/home/tools/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar

RUN cd /home && wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
RUN cd /home/tools/strelka-2.9.2.centos6_x86_64 #&& tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2
RUN cd /home/tools/strelka-2.9.2.centos6_x86_64 && bash bin/runStrelkaSomaticWorkflowDemo.bash
RUN cd /home/tools/strelka-2.9.2.centos6_x86_64 && bash bin/runStrelkaGermlineWorkflowDemo.bash

RUN cd /home/tools && tar xvf R-3.6.3.tar.gz
RUN cd /home/tools/R-3.6.3 && ./configure
RUN cd /home/tools/R-3.6.3 && make

RUN cd /home/tools && wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
RUN cd /home/tools && tar -xf bcftools-1.10.2.tar.bz2
RUN cd /home/tools/bcftools-1.10.2 && ./configure
RUN cd /home/tools/bcftools-1.10.2 && make
RUN cd /home/tools/bcftools-1.10.2 && make install

RUN apt-get update
RUN apt-get install -y default-jre
RUN apt-get install -y default-jdk
RUN echo 'JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java"' >> /etc/environment



