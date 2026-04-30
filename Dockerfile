# Builds the thermoembo C++ FEM solver and bundles the Python pipeline
# (run_thermoembo.py) so the entire workflow runs inside one container.
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
# Run the full pipeline (mount a directory containing the labeled NIfTI):
#   docker run --rm -v $(pwd):/data thermoembo \
#     /data/PreTxArtLoRes.vessellabel.nii.gz \
#     --out-dir /data/thermoembo_run --steps 5
#
# Prepare meshes only (skip simulation):
#   docker run --rm -v $(pwd):/data thermoembo \
#     /data/label.nii.gz --out-dir /data/out --mesh-only
#
# Drop to a shell (override entrypoint):
#   docker run --rm -it --entrypoint /bin/bash -v $(pwd):/data thermoembo

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
#
# Debug / IDE tools:
# gdb / cgdb:            debugger and curses TUI frontend
# valgrind:              memory error detector and profiler
# vim / exuberant-ctags / cscope: editor + code-navigation index generators
# tmux:                  terminal multiplexer (split panes for gdb + editor)
# strace / ltrace:       syscall / library-call tracing
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ gfortran \
        make wget curl python \
        ca-certificates \
        cmake m4 file \
        openmpi-bin libopenmpi-dev \
        libvtk5-dev \
        liblapack-dev libblas-dev \
        zlib1g-dev \
        gdb cgdb valgrind \
        vim exuberant-ctags cscope \
        tmux strace ltrace \
        less man-db \
    && rm -rf /var/lib/apt/lists/*

# ── Vim configuration ────────────────────────────────────────────────────────
RUN cat > /root/.vimrc <<'VIMRC'
set number ruler showcmd laststatus=2
set expandtab tabstop=4 shiftwidth=4 autoindent
set hlsearch incsearch ignorecase smartcase
set mouse=a
set wildmenu
set background=dark
syntax on

" ctags — search from current dir up to /work/tutorials
set tags=./tags;/work/tutorials

" cscope — auto-connect if database exists beside the source
if has("cscope")
  set csprg=/usr/bin/cscope
  set csto=0 cst nocsverb
  if filereadable("cscope.out")
    cs add cscope.out
  elseif filereadable("/work/tutorials/cscope.out")
    cs add /work/tutorials/cscope.out
  endif
  set csverb
  nmap <C-\>s :cs find s <C-R>=expand("<cword>")<CR><CR>
  nmap <C-\>g :cs find g <C-R>=expand("<cword>")<CR><CR>
  nmap <C-\>c :cs find c <C-R>=expand("<cword>")<CR><CR>
  nmap <C-\>f :cs find f <C-R>=expand("<cfile>")<CR><CR>
endif

" File tree with :Explore / :e .
let g:netrw_liststyle=3 g:netrw_banner=0
VIMRC

# ── PETSc 3.10.2 ─────────────────────────────────────────────────────────────
ENV PETSC_VERSION=3.10.2
ENV PETSC_DIR=/opt/petsc/petsc-3.10.2
ENV PETSC_ARCH=3.10.2-xenial-gcc-5.4.0-opt

RUN mkdir -p /opt/petsc \
 && wget -q \
        https://github.com/petsc/petsc/archive/refs/tags/v${PETSC_VERSION}.tar.gz \
        -O /tmp/petsc.tar.gz \
 && tar -xzf /tmp/petsc.tar.gz -C /opt/petsc \
 && rm /tmp/petsc.tar.gz

# NetCDF 4.5.0 — PETSc 3.10.2 hardcodes a dead Unidata FTP URL.
# The GitHub archive extracts to netcdf-c-4.5.0/ but PETSc infers the
# directory name from the tarball (netcdf-4.5.0), so we repack it.
RUN curl -sL https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.5.0.tar.gz \
    | tar xzf - -C /tmp \
 && mv /tmp/netcdf-c-4.5.0 /tmp/netcdf-4.5.0 \
 && tar czf /tmp/netcdf-4.5.0.tar.gz -C /tmp netcdf-4.5.0 \
 && rm -rf /tmp/netcdf-4.5.0

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
        --with-zlib-dir=/usr \
        --download-hdf5=1 \
        --download-pnetcdf=1 \
        --download-netcdf=/tmp/netcdf-4.5.0.tar.gz \
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

RUN cd /work/tutorials && make thermoembo thermoembo1d thermoembo-debug thermoembo1d-debug

# ── Python environment (Miniforge) ───────────────────────────────────────────
# Miniforge defaults to conda-forge and ships mamba for faster solving.
# vtk is capped at 9.2.* — the last series whose conda-forge builds target
# the cos7 sysroot (GLIBC 2.17).  VTK 9.3+ migrated to cos8 (GLIBC 2.28).
# pyvista <0.40 is the last branch compatible with vtk 9.2.
RUN wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
        -O /tmp/miniforge.sh \
 && bash /tmp/miniforge.sh -b -p /opt/conda \
 && rm /tmp/miniforge.sh
ENV PATH=/opt/conda/bin:${PATH}

# conda handles lightweight packages; vtk/pyvista/tetgen installed via pip
# because the vtk conda-forge package pulls in a huge dependency tree that
# causes the libmamba SAT solver to stall on the python=3.11 constraint.
# pip ships a self-contained vtk wheel (manylinux2014, glibc ≥2.17) that
# installs in seconds and works on Ubuntu 16.04 (glibc 2.23).
# tetgen capped at 0.6.4 and pymeshfix capped at 0.16.3 — the last releases
# that ship manylinux2014 (glibc ≥2.17) pre-built wheels.  Newer versions
# (tetgen ≥0.7, pymeshfix ≥0.17) switched to nanobind (C++17) and publish
# manylinux_2_28 wheels only, which require glibc ≥2.28 (not available on
# Ubuntu 16.04 / glibc 2.23).
RUN mamba install -y python=3.10 \
        nibabel scipy scikit-image meshio netcdf4 "numpy<2" pandas \
 && pip install "vtk==9.2.6" "pyvista==0.38.6" "pymeshfix==0.16.3" "tetgen==0.6.4" \
 && conda clean -afy

# ── Pipeline scripts ──────────────────────────────────────────────────────────
COPY pipeline/ /work/pipeline/
ENV PYTHONPATH=/work/pipeline:${PYTHONPATH}

# ── Runtime environment ───────────────────────────────────────────────────────
# Allow mpirun as root inside the container (common in HPC containers)
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /data
ENTRYPOINT ["python", "/work/pipeline/run_thermoembo.py"]
CMD ["--help"]
