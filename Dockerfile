# DOCKER-VERSION 1.12.3

FROM genomehubs/easy-import:89
MAINTAINER Richard Challis / Lepbase contact@lepbase.org

ENV TERM xterm
ENV DEBIAN_FRONTEND noninteractive

USER root
RUN cpanm Parallel::ForkManager

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y libbz2-dev

WORKDIR /
RUN git clone git://github.com/samtools/htslib.git
RUN git clone git://github.com/samtools/bcftools.git
WORKDIR /bcftools
RUN make

USER eguser
WORKDIR /ensembl
#RUN git clone -b release/89 https://github.com/ensembl/ensembl-tools

RUN git clone https://github.com/adamsardar/perl-libs-custom.git

WORKDIR /ensembl/easy-import
ARG  cachebuster=0b7ad45c8
RUN  git pull origin develop && git submodule update --recursive

COPY startup.sh /import/

ENV PATH $PATH:/bcftools
ENV PERL5LIB $PERL5LIB:/ensembl/perl-libs-custom/EnsemblAPI/ensembl-variation/scripts/import
ENV HOME tmp

CMD /import/startup.sh $FLAGS
