
all: thermoembo
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

thermoembo: thermoembo.o chkopts
	-${CLINKER} -o thermoembo-${PETSC_ARCH} thermoembo.o ${PETSC_TS_LIB} ${VTK_LIB}

phasemodel: phasemodel.o chkopts
	-${CLINKER} -o phasemodel phasemodel.o ${PETSC_TS_LIB} ${VTK_LIB}

exac: exac.o chkopts
	-${CLINKER} -o exac exac.o ${PETSC_TS_LIB} ${VTK_LIB}

tags:
	ctags -R --languages=c,c++ --langmap=c++:+.txx --langmap=c++:+.cu --langmap=c++:+.cl $(PETSC_DIR)  /usr/include/vtk-5.10/ .
