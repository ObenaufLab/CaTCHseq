FROM rust:1-slim-bookworm AS builder

WORKDIR /src
ADD source.tar.gz /src
COPY Cargo.lock /src/Cargo.lock
RUN cargo build --release --locked

FROM python:3.12-slim-bookworm

RUN pip install --no-cache-dir cutadapt==5.0 PyYAML==6.0.2

COPY --from=builder /src/target/release/pipspeak /usr/local/bin/pipspeak

RUN pipspeak --version \
    && cutadapt --version \
    && python3 -c "import yaml; assert yaml.__version__ == '6.0.2', yaml.__version__"

LABEL maintainer="joerg.fallmann@imp.ac.at"
LABEL version="e2524c7"
LABEL description="pipspeak PIPseq to 10X-compatible FASTQ preprocessing"
