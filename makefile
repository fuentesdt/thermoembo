SHELL := /bin/bash
C3DEXE=c3d
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
	$(C3DEXE) -verbose $< -smooth 1x1x0vox -hessobj 1  .5  3   -o $@

$(WORKDIR)/vessel.nii.gz: $(WORKDIR)/HESSOBJ14.nii.gz
	$(C3DEXE) $< -thresh 2 inf 1 0 $(WORKDIR)/mask.nii.gz -thresh 1 1 1 0 -multiply -dilate 1 3x3x3vox  -erode 1 3x3x3vox  -o $@

$(WORKDIR)/hesstmp.nii.gz: 
	$(C3DEXE)  $< -thresh 10 inf 1 0 $(WORKDIR)/mask.nii.gz -thresh 1 1 1 0 -multiply -o $(WORKDIR)/hesstmp.nii.gz
	$(C3DEXE) $(WORKDIR)/anatomy.nii.gz $(WORKDIR)/hesstmp.nii.gz -lstat

$(WORKDIR)/anatomyblur.nii.gz:
	$(C3DEXE) -verbose $(WORKDIR)/anatomy.nii.gz -smooth 1x1x0vox   -o $@

$(WORKDIR)/vesseltmp.nii.gz:
	$(C3DEXE) -verbose $(WORKDIR)/anatomyblur.nii.gz -thresh 120 inf 1 0 $(WORKDIR)/mask.nii.gz -thresh 1 1 1 0 -multiply -comp -thresh 1 15 1 0 -o $@

$(WORKDIR)/vesselnew.nii.gz: $(WORKDIR)/vesseltmp.nii.gz
	$(C3DEXE) -verbose $(WORKDIR)/anatomy.nii.gz  -thresh 140 inf 1 -1  $(WORKDIR)/hesstmp.nii.gz -replace 0 1 1 -1 -levelset-curvature 1.0 -levelset 10  -thresh -inf 0 1 0 -o $@

datalocation/radiomicsout.csv: datalocation/radiomics.csv
	python /opt/apps/pyradiomics/radiomics/scripts/commandlinebatch.py $< $@ -p Params.yaml  -l 1  -l 2 -v  5

viewseg:
	vglrun /opt/apps/itksnap/itksnap-3.6.0-20170401-Linux-x86_64/bin/itksnap -g datalocation/anatomyblur.nii.gz -s datalocation/vesseltmp.nii.gz  -o  datalocation/HESSOBJ.nii.gz  datalocation/hesstmp.nii.gz

viewhess:
	vglrun /opt/apps/itksnap/itksnap-3.6.0-20170401-Linux-x86_64/bin/itksnap -g datalocation/anatomy.nii.gz -s datalocation/hesstmp.nii.gz  -o datalocation/HESSOBJ.nii.gz
view:
	vglrun /opt/apps/itksnap/itksnap-3.6.0-20170401-Linux-x86_64/bin/itksnap -g datalocation/anatomy.nii.gz -s datalocation/label.nii.gz -o datalocation/HESSOBJ*.nii.gz
