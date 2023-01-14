# syntax=docker/dockerfile:1
FROM python:3.9-buster
RUN python --version
WORKDIR /pipeline
#RUN python -m pip install --upgrade pip
COPY req.txt req.txt
RUN pip3 install -r req.txt
COPY . . 
CMD /bin/bash
