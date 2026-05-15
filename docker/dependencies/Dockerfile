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
      gfortran \
      ca-certificates \
      clang \
      && \
    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb && \
    dpkg -i cuda-keyring_1.1-1_all.deb && \
    rm cuda-keyring_1.1-1_all.deb && \
    apt-get update && apt-get install -y --no-install-recommends \
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
      libvtk9-dev \
      cuda-nvcc-12-8 \
      cuda-cudart-dev-12-8 \
      cuda-driver-dev-12-8 \
      libcublas-dev-12-8 \
      libcusolver-dev-12-8 \
      libcusparse-dev-12-8 \
      libcurand-dev-12-8 \
      libnvjitlink-12-8 \
      libnvjitlink-dev-12-8 \
      cuda-nvtx-12-8 \
  && rm -rf /var/lib/apt/lists/*

# Create directory for dependencies
ARG NPROCS=12
ENV NPROCS=$NPROCS \
    DEPS_ROOT="/opt/mirco-dependencies"
RUN mkdir -p ${DEPS_ROOT}

COPY dependencies /dependencies

# Make `nvcc` available as an environment variable
ENV PATH=/usr/local/cuda/bin:${PATH}

# Install Kokkos and Kokkos-Kernels
RUN /dependencies/kokkos/install.sh ${DEPS_ROOT}
RUN /dependencies/kokkos-kernels/install.sh ${DEPS_ROOT}

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
