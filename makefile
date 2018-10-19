SHELL := /bin/bash
C3DEXE=/rsrch2/ip/dtfuentes/bin/c3d
WORKDIR=datalocation


$(WORKDIR)/HESSOBJ17.nii.gz: $(WORKDIR)/anatomy.nii.gz
	for idvar in `seq 1 1`; do $(C3DEXE) -verbose $<  -hessobj $$idvar .5 .5   -o $(WORKDIR)/HESSOBJ$${idvar}0.nii.gz; done
	for idvar in `seq 1 1`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  1  1   -o $(WORKDIR)/HESSOBJ$${idvar}1.nii.gz; done
	for idvar in `seq 1 1`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  2  2   -o $(WORKDIR)/HESSOBJ$${idvar}2.nii.gz; done
	for idvar in `seq 1 1`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  3  3   -o $(WORKDIR)/HESSOBJ$${idvar}3.nii.gz; done
	for idvar in `seq 1 1`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  4  4   -o $(WORKDIR)/HESSOBJ$${idvar}4.nii.gz; done
	for idvar in `seq 1 1`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  5  5   -o $(WORKDIR)/HESSOBJ$${idvar}5.nii.gz; done
	for idvar in `seq 1 1`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  6  6   -o $(WORKDIR)/HESSOBJ$${idvar}6.nii.gz; done
	for idvar in `seq 1 1`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  7  7   -o $(WORKDIR)/HESSOBJ$${idvar}7.nii.gz; done
$(WORKDIR)/vesselcomp.nii.gz: $(WORKDIR)/HESSOBJ17.nii.gz
	$(C3DEXE) $< -thresh 1 inf 1 0 -comp  -o $@
