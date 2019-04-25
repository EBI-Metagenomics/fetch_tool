FROM python:3.6-alpine
LABEL maintainer="Miguel Boland <mdb@ebi.ac.uk>"

# Python build dependencies

COPY src /src
COPY setup.py /
COPY requirements.txt /
COPY __init__.py /

RUN apk --no-cache add --virtual .builddeps gcc gfortran musl-dev

RUN pip install -r requirements.txt
RUN pip --no-cache-dir -v install .

RUN apk del .builddeps && rm -rf /root/.cache