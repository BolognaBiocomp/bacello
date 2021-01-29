# Base Image
FROM python:3.8-slim-buster

# Metadata
LABEL base.image="python:3.8-slim-buster"
LABEL version="1.0"
LABEL software="BaCelLo"
LABEL software.version="202001"
LABEL description="an open source software tool to predict subcellular localzation of eukaryotic proteins"
LABEL website="https://github.com/BolognaBiocomp/bacello"
LABEL documentation="https://github.com/BolognaBiocomp/bacello"
LABEL license="GNU GENERAL PUBLIC LICENSE Version 3"
LABEL tags="Proteomics"
LABEL maintainer="Castrense Savojardo <castrense.savojardo2@unibo.it>"

WORKDIR /usr/src/bacello

COPY requirements.txt .

RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    useradd -m bacello

USER bacello

COPY . .

WORKDIR /data/

ENV BACELLO_HOME='/usr/src/bacello' PATH=/usr/src/bacello:$PATH

ENTRYPOINT ["/usr/src/bacello/bacello.py"]
