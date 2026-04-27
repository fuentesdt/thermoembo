# Builds and runs tutorials/thermoembo and tutorials/thermoembo1d
#
# Version pins match the repository's PETSC_ARCH string:
#   3.10.2-xenial-gcc-5.4.0-opt
#   => PETSc 3.10.2, Ubuntu 16.04 (Xenial), GCC 5.4 (default on Xenial), optimised
#
# VTK 5.10 is the only Ubuntu-packaged VTK version that ships the exact
# /usr/lib/libvtk*.so.5.10 paths hard-coded in tutorials/makefile.
#
# Build:
#   docker build -t thermoembo .
#
# Run (mount a working directory containing mesh + image files):
#   docker run --rm -v $(pwd):/data thermoembo \
#     mpirun --allow-run-as-root -n 4 \
#       /work/tutorials/thermoembo-3.10.2-xenial-gcc-5.4.0-opt \
#       -dim 3 -mesh /data/mytetmesh.2.exo -vtk /data/output.vtk \
#       -temp_petscspace_degree 1 -pres_petscspace_degree 1 \
#       -ts_type beuler -ts_max_steps 100 -ts_dt 1.e-1

FROM ubuntu:16.04

# Suppress interactive prompts during apt
ENV DEBIAN_FRONTEND=noninteractive

# ── System packages ──────────────────────────────────────────────────────────
# gcc/g++/gfortran: GCC 5.4.0 (Ubuntu 16.04 default)
# libvtk5-dev:      VTK 5.10 headers + /usr/lib/libvtk*.so.5.10 (exact paths
#                   expected by tutorials/makefile)
# libopenmpi-dev:   MPI headers and mpicc/mpicxx/mpif90 wrappers
# liblapack-dev / libblas-dev: BLAS+LAPACK for PETSc linear algebra
# zlib1g-dev:       required by NetCDF (downloaded by PETSc)
# python:           Python 2.7 — required by PETSc's configure script
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ gfortran \
        make wget python \
        openmpi-bin libopenmpi-dev \
        libvtk5-dev \
        liblapack-dev libblas-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ── PETSc 3.10.2 ─────────────────────────────────────────────────────────────
ENV PETSC_VERSION=3.10.2
ENV PETSC_DIR=/opt/petsc/petsc-3.10.2
ENV PETSC_ARCH=3.10.2-xenial-gcc-5.4.0-opt

RUN mkdir -p /opt/petsc \
 && wget -q \
        https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VERSION}.tar.gz \
        -O /tmp/petsc.tar.gz \
 && tar -xzf /tmp/petsc.tar.gz -C /opt/petsc \
 && rm /tmp/petsc.tar.gz

# Configure options:
#   --download-netcdf / --download-exodusii: required for DMPlexCreateFromFile
#     to load .exo (Exodus II) meshes; no Ubuntu package for ExodusII exists
#   --with-blaslapack-lib: explicit because Ubuntu 16.04 BLAS/LAPACK live in
#     /usr/lib/x86_64-linux-gnu, not /usr/lib, so the default dir search fails
#   --with-shared-libraries=0: static PETSc avoids LD_LIBRARY_PATH headaches
#     and matches how the tutorials makefile links (${PETSC_TS_LIB} is static)
RUN cd ${PETSC_DIR} && ./configure \
        PETSC_ARCH=${PETSC_ARCH} \
        --with-cc=mpicc \
        --with-cxx=mpicxx \
        --with-fc=mpif90 \
        --with-debugging=0 \
        --COPTFLAGS='-O2' \
        --CXXOPTFLAGS='-O2' \
        --FOPTFLAGS='-O2' \
        --with-blaslapack-lib='-L/usr/lib/x86_64-linux-gnu -llapack -lblas -lgfortran -lm' \
        --download-netcdf=1 \
        --download-exodusii=1 \
        --with-shared-libraries=0 \
 && make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all

# ── Source ────────────────────────────────────────────────────────────────────
WORKDIR /work
COPY tutorials/ /work/tutorials/

# ── Build ─────────────────────────────────────────────────────────────────────
# VTK 5.10 runtime libs are in /usr/lib (Ubuntu 16.04 package layout).
# OpenMPI runtime libs need LD_LIBRARY_PATH for mpirun at link time.
ENV LD_LIBRARY_PATH=/usr/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:${LD_LIBRARY_PATH}

RUN cd /work/tutorials && make thermoembo thermoembo1d

# ── Runtime environment ───────────────────────────────────────────────────────
# Allow mpirun as root inside the container (common in HPC containers)
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /data
CMD ["/bin/bash"]
