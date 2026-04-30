
all: thermoembo thermoembo1d thermoembo-debug thermoembo1d-debug
CFLAGS          =
FFLAGS          =
CPPFLAGS	= -I/usr/include/vtk-5.10
FPPFLAGS        =
LOCDIR          = .
EXAMPLESC       = exac.c ex45.c
EXAMPLESF       = ex1f.F ex22f.F ex22f_mf.F90
MANSEC          = TS
DIRS            = phasefield advection-diffusion-reaction eimex power_grid network
CLEANFILES      = SA-data/*
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

VTK_LIB = -rdynamic /usr/lib/libvtkVolumeRendering.so.5.10 /usr/lib/libvtkIO.so.5.10 /usr/lib/libvtkRendering.so.5.10 /usr/lib/libvtkIO.so.5.10 /usr/lib/libvtkGraphics.so.5.10 /usr/lib/libvtkImaging.so.5.10 /usr/lib/libvtkFiltering.so.5.10 /usr/lib/libvtkCommon.so.5.10 /usr/lib/libvtksys.so.5.10 /usr/lib/libvtkHybrid.so.5.10

# thermoembo*.c include C++ VTK headers; override PETSc's default .c.o rule
# (which uses mpicc) with an explicit C++ rule using ${CXX} = mpicxx.
thermoembo.o: thermoembo.c
	${CXX} ${CXXFLAGS} ${PETSC_CCPPFLAGS} ${CPPFLAGS} -c -o $@ $<

thermoembo1d.o: thermoembo1d.c
	${CXX} ${CXXFLAGS} ${PETSC_CCPPFLAGS} ${CPPFLAGS} -c -o $@ $<

thermoembo: thermoembo.o chkopts
	-${CLINKER} -o thermoembo-${PETSC_ARCH} thermoembo.o ${PETSC_TS_LIB} ${VTK_LIB}

thermoembo1d: thermoembo1d.o chkopts
	-${CLINKER} -o thermoembo1d-${PETSC_ARCH} thermoembo1d.o ${PETSC_TS_LIB} ${VTK_LIB}

# ── Debug builds (-O0 -g, no inlining) ───────────────────────────────────────
# CXXFLAGS is omitted so PETSc's -O2 optimisation flags do not override -O0.
# PETSC_CCPPFLAGS provides all necessary include paths and defines.
DEBUG_FLAGS = -O0 -g3 -fno-omit-frame-pointer

thermoembo-debug.o: thermoembo.c
	${CXX} ${DEBUG_FLAGS} ${PETSC_CCPPFLAGS} ${CPPFLAGS} -c -o $@ $<

thermoembo1d-debug.o: thermoembo1d.c
	${CXX} ${DEBUG_FLAGS} ${PETSC_CCPPFLAGS} ${CPPFLAGS} -c -o $@ $<

thermoembo-debug: thermoembo-debug.o chkopts
	-${CLINKER} -o thermoembo-${PETSC_ARCH}-debug thermoembo-debug.o ${PETSC_TS_LIB} ${VTK_LIB}

thermoembo1d-debug: thermoembo1d-debug.o chkopts
	-${CLINKER} -o thermoembo1d-${PETSC_ARCH}-debug thermoembo1d-debug.o ${PETSC_TS_LIB} ${VTK_LIB}

phasemodel: phasemodel.o chkopts
	-${CLINKER} -o phasemodel phasemodel.o ${PETSC_TS_LIB} ${VTK_LIB}

exac: exac.o chkopts
	-${CLINKER} -o exac exac.o ${PETSC_TS_LIB} ${VTK_LIB}

tags:
	ctags -R --languages=c,c++ --langmap=c++:+.txx --langmap=c++:+.cu --langmap=c++:+.cl $(PETSC_DIR)  /usr/include/vtk-5.10/ .
