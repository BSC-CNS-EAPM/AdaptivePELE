FROM python:3.7-slim as builder
LABEL maintainer="cescgina@gmail.com"

RUN apt-get update \
    && apt-get install gcc -y \
    && apt-get clean

WORKDIR /AdaptivePELE

COPY requirements.txt /AdaptivePELE

RUN pip install --no-cache-dir --requirement requirements.txt

COPY . /AdaptivePELE

RUN python setup.py install

FROM python:3.7-slim as adaptivepele

COPY --from=builder /usr/local/lib/python3.7/site-packages /usr/local/lib/python3.7/site-packages
