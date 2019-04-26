FROM python:3.6
LABEL maintainer="Miguel Boland <mdb@ebi.ac.uk>"

# Python build dependencies

COPY src /src
COPY setup.py /
COPY requirements.txt /
COPY __init__.py /

RUN pip install -r requirements.txt
RUN pip --no-cache-dir -v install .

