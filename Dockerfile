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

COPY coherens/ /opt/coherens/
COPY setups/ /opt/coherens/setups/
COPY docker/ /docker/

ENV COHERENS_ROOT=/opt/coherens \
    COHERENS_FLAGS_FILE=/docker/coherensflags_docker.cmp

WORKDIR /workspace
RUN mkdir -p /workspace/runs

ENTRYPOINT ["/docker/run_setup.sh"]
CMD ["seamount_smoke"]
