FROM python:3.7-slim
LABEL maintainer="cescgina@gmail.com"

RUN apt-get update \
    && apt-get install gcc -y \
    && apt-get clean

WORKDIR /AdaptivePELE

COPY requirements.txt /AdaptivePELE

RUN pip install --no-cache-dir --requirement requirements.txt

COPY . /AdaptivePELE

RUN python setup.py build_ext --inplace && rm -r build
