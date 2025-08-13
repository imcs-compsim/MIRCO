ARG BASE_IMAGE=ubuntu:24.04
FROM ${BASE_IMAGE}
LABEL org.opencontainers.image.description="Image containing all the dependencies required for building and testing MIRCO"
LABEL org.mirco.project=MIRCO

# Prevents tzdata asking for user feedback
ENV DEBIAN_FRONTEND=noninteractive

USER root

# Set locale information:. region and timezone
RUN apt-get update && apt-get install -y locales && \
    localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8 && \
    locale-gen en_US.UTF-8 && \
    rm -rf /var/lib/apt/lists/*
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

RUN apt-get update && apt-get install -y \
      build-essential \
      cmake \
      git \
      libglu1-mesa \
      python3 \
      sudo \
      unzip \
      vim \
      wget \
      && \
    apt-get update && apt-get install -y \
      libblas-dev \
      libboost-all-dev \
      libcln-dev \
      libhdf5-dev \
      libhdf5-openmpi-dev \
      libnetcdf-dev \
      libfftw3-dev \
      lld \
      python3 \
      python3-venv \
      python-is-python3 \
      liblapack-dev \
      ninja-build \
      libomp-dev \
  && rm -rf /var/lib/apt/lists/*

# Create directory for dependencies
ARG NPROCS=12
ENV NPROCS=$NPROCS \
    INSTALL_DIR="/opt/mirco-dependencies/"
RUN mkdir -p ${INSTALL_DIR}

COPY dependencies /dependencies

# Install Trilinos 2022_Q1
RUN /dependencies/trilinos/install.sh ${INSTALL_DIR}/trilinos

# Install Kokkos and Kokkos-Kernels
RUN /dependencies/kokkos/install.sh ${INSTALL_DIR}/kokkos
RUN /dependencies/kokkos-kernels/install.sh ${INSTALL_DIR}/kokkos-kernels

# Add dependencies hash
# The label is added at the end because the label causes a cache miss
ARG DEPENDENCIES_HASH
LABEL org.mirco.dependencies_hash="${DEPENDENCIES_HASH}"
ENV DEPENDENCIES_HASH=${DEPENDENCIES_HASH}

# add and enable the default user
ENV USER=user
RUN adduser --disabled-password --shell '/usr/bin/bash' --gecos '' $USER
RUN adduser $USER sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

#make sure everything is in place
RUN chown -R $USER:$USER /home/$USER
USER $USER
ENV HOME=/home/$USER
ENV USER=$USER
# see https://github.com/open-mpi/ompi/issues/4948
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
WORKDIR $HOME
