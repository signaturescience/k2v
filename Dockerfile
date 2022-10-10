FROM alpine:3.15

ARG BCFTOOLS_VERSION="1.9"

ADD src /src
WORKDIR /src
RUN apk update && apk add --no-cache gcc g++ make libc-dev ncurses-dev zlib-dev xz-dev bzip2-dev R
RUN tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    make -j4 && \
    make install && \
    cd /src && \
    rm -rf bcftools-${BCFTOOLS_VERSION} bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN tar -xjf htslib-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd htslib-${BCFTOOLS_VERSION} && \
    make -j4 && \
    make install && \
    cd /src && \
    rm -rf htslib-${BCFTOOLS_VERSION} htslib-${BCFTOOLS_VERSION}.tar.bz2

ADD assets /assets
ADD testdata /testdata
ADD scripts /scripts

WORKDIR /data

ENTRYPOINT ["/scripts/k2v.sh"]
