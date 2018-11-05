SHELL := /bin/bash
C3DEXE=/rsrch2/ip/dtfuentes/bin/c3d
WORKDIR=datalocation


$(WORKDIR)/HESSOBJ27.nii.gz: $(WORKDIR)/anatomy.nii.gz
	for idvar in `seq 0 2`; do $(C3DEXE) -verbose $<  -hessobj $$idvar .5 .5   -o $(WORKDIR)/HESSOBJ$${idvar}0.nii.gz; done
	for idvar in `seq 0 2`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  1  1   -o $(WORKDIR)/HESSOBJ$${idvar}1.nii.gz; done
	for idvar in `seq 0 2`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  2  2   -o $(WORKDIR)/HESSOBJ$${idvar}2.nii.gz; done
	for idvar in `seq 0 2`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  3  3   -o $(WORKDIR)/HESSOBJ$${idvar}3.nii.gz; done
	for idvar in `seq 0 2`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  4  4   -o $(WORKDIR)/HESSOBJ$${idvar}4.nii.gz; done
	for idvar in `seq 0 2`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  5  5   -o $(WORKDIR)/HESSOBJ$${idvar}5.nii.gz; done
	for idvar in `seq 0 2`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  6  6   -o $(WORKDIR)/HESSOBJ$${idvar}6.nii.gz; done
	for idvar in `seq 0 2`; do $(C3DEXE) -verbose $<  -hessobj $$idvar  7  7   -o $(WORKDIR)/HESSOBJ$${idvar}7.nii.gz; done

$(WORKDIR)/HESSOBJ.nii.gz: $(WORKDIR)/anatomy.nii.gz
	$(C3DEXE) -verbose $<  -hessobj 1  1  5   -o $@

$(WORKDIR)/vessel.nii.gz: $(WORKDIR)/HESSOBJ14.nii.gz
	$(C3DEXE) $< -thresh 2 inf 1 0 $(WORKDIR)/mask.nii.gz -thresh 1 1 1 0 -multiply -dilate 1 3x3x3vox  -erode 1 3x3x3vox  -o $@

datalocation/radiomicsout.csv: datalocation/radiomics.csv
	python /opt/apps/pyradiomics/radiomics/scripts/commandlinebatch.py $< $@ -p Params.yaml  -l 1  -l 2 -v  5

viewseg:
	vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g datalocation/anatomy.nii.gz -s datalocation/vessel.nii.gz  -o datalocation/HESSOBJ14.nii.gz
viewhess:
	vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g datalocation/anatomy.nii.gz -s datalocation/mask.nii.gz  -o datalocation/HESSOBJ.nii.gz
view:
	vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g datalocation/anatomy.nii.gz -s datalocation/label.nii.gz -o datalocation/HESSOBJ*.nii.gz
