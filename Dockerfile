FROM centos:centos7
LABEL maintainer="Miguel Boland <mdb@ebi.ac.uk>"

# Python build dependencies

COPY src /src
COPY setup.py /
COPY requirements.txt /
COPY __init__.py /

RUN yum -y update; yum clean all

RUN yum -y install epel-release gcc bzip2 git wget; yum clean all

RUN yum -y install python36 python36-devel python36-setuptools
RUN easy_install-3.6 pip
RUN pip3 install -U pip virtualenv setuptools

RUN pip3 install -r requirements.txt
RUN python3 setup.py install

