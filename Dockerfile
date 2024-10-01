# Use python 3.12
FROM python:3.12-slim

ENV PYTHONHASHSEED=0
ENV PYTHONUNBUFFERED=1

# turn off built-in Python multithreading
ENV MKL_NUM_THREADS=1
ENV NUMEXPR_NUM_THREADS=1
ENV OMP_NUM_THREADS=1

# setup the working directory
RUN mkdir /code && \
    mkdir /code/RVData && \
    mkdir /data && \
    mkdir /outputs && \
    apt-get --yes update && \
    apt install build-essential -y --no-install-recommends && \
    apt-get install --yes git vim emacs nano && \
    /usr/local/bin/python -m pip install --upgrade pip && \
    cd /code/RVData && \
    mkdir -p logs && \
	mkdir -p outputs

# Set the working directory to RVData
WORKDIR /code/RVData
RUN git config --global --add safe.directory /code/RVData

ADD requirements.txt /code/RVData/
RUN pip3 install -r requirements.txt
