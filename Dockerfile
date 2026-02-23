FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       gfortran \
       git \
       libnetcdf-dev \
       libnetcdff-dev \
       libopenmpi-dev \
       make \
       openmpi-bin \
       time \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/2026_coherens_container
COPY . .

ENV COHERENS_ROOT=/opt/2026_coherens_container/coherens \
    COHERENS_FLAGS_FILE=/opt/2026_coherens_container/docker/coherensflags_docker.cmp

WORKDIR /workspace
RUN mkdir -p /workspace/runs

ENTRYPOINT ["/opt/2026_coherens_container/docker/run_setup.sh"]
CMD ["tutorial/river"]
