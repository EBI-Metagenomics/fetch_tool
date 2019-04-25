FROM alpine:3.9
RUN echo "http://dl-8.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories
RUN apk --no-cache --update-cache add gcc gfortran python python-dev py-pip build-base wget freetype-dev libpng-dev openblas-dev
RUN ln -s /usr/include/locale.h /usr/include/xlocale.h

LABEL maintainer="Miguel Boland <mdb@ebi.ac.uk>"

# Python build dependencies

COPY src /src
COPY setup.py /
COPY requirements.txt /
COPY __init__.py /

RUN pip install .
