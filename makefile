SHELL := /bin/bash
C3DEXE=c3d
WORKDIR=datalocation

UIDLIST=ZZ17S003 ZZ17S004 ZZ17S005 ZZ17S006 ZZ17S007

# keep tmp files
.SECONDARY: 

hess:      $(addprefix $(WORKDIR)/,$(addsuffix /hessobj.nii.gz,$(UIDLIST)))  
hessseg:   $(addprefix $(WORKDIR)/,$(addsuffix /hessseg.nii.gz,$(UIDLIST)))  
vessel:    $(addprefix $(WORKDIR)/,$(addsuffix /vessel.nii.gz,$(UIDLIST)))  
mask:      $(addprefix $(WORKDIR)/,$(addsuffix /mask,$(UIDLIST)))  
viewhess:  $(addprefix $(WORKDIR)/,$(addsuffix /viewhess,$(UIDLIST)))  
viewvess:  $(addprefix $(WORKDIR)/,$(addsuffix /viewvess,$(UIDLIST)))  

# https://sourceforge.net/p/c3d/git/ci/master/tree/adapters/HessianObjectness.cxx#l66
$(WORKDIR)/%/hessobj.nii.gz: $(WORKDIR)/%/anatomy.nii.gz
	$(C3DEXE) -verbose $<  -hessobj 1  .1  3   -o $@

$(WORKDIR)/%/mask: $(WORKDIR)/%/anatomy.nii.gz
	if [ ! -f $@.nii.gz  ] ; then c3d $<  -scale 0 -type uchar -o $@.nii.gz;  fi
	vglrun /opt/apps/itksnap/itksnap-3.6.0-20170401-Linux-x86_64/bin/itksnap -g $< -s $@.nii.gz -o $(WORKDIR)/$*/hessobj.nii.gz

$(WORKDIR)/%/hessseg.nii.gz: $(WORKDIR)/%/hessobj.nii.gz  $(WORKDIR)/%/mask.nii.gz
	$(C3DEXE)  $< -thresh 25 inf 1 0 $(word 2,$^) -thresh 1 1 1 0 -multiply -o $@
	$(C3DEXE) $(WORKDIR)/$*/anatomy.nii.gz $@ -lstat

$(WORKDIR)/%/vessel.nii.gz:
	$(C3DEXE) -verbose $(WORKDIR)/$*/anatomy.nii.gz -thresh 120 inf 1 0 $(WORKDIR)/$*/mask.nii.gz -thresh 1 1 1 0 -multiply -comp -thresh 1 15 1 0 -o $@

$(WORKDIR)/%/viewpre: $(WORKDIR)/%/preone.nii.gz   $(WORKDIR)/%/pretwo.nii.gz
	vglrun /opt/apps/itksnap/itksnap-3.6.0-20170401-Linux-x86_64/bin/itksnap -g  $< -s $(WORKDIR)/$*/preseg.nii.gz  -o $(word 2,$^) 

$(WORKDIR)/%/viewhess: $(WORKDIR)/%/anatomy.nii.gz
	vglrun /opt/apps/itksnap/itksnap-3.6.0-20170401-Linux-x86_64/bin/itksnap -g  $< -s $(WORKDIR)/$*/hessseg.nii.gz  -o $(WORKDIR)/$*/hessobj.nii.gz

$(WORKDIR)/%/viewvess: $(WORKDIR)/%/anatomy.nii.gz
	vglrun /opt/apps/itksnap/itksnap-3.6.0-20170401-Linux-x86_64/bin/itksnap -g  $< -s $(WORKDIR)/$*/vessel.nii.gz  -o $(WORKDIR)/$*/hessobj.nii.gz


$(WORKDIR)/vessel.nii.gz: $(WORKDIR)/HESSOBJ14.nii.gz
	$(C3DEXE) $< -thresh 2 inf 1 0 $(WORKDIR)/mask.nii.gz -thresh 1 1 1 0 -multiply -dilate 1 3x3x3vox  -erode 1 3x3x3vox  -o $@

$(WORKDIR)/vesselnew.nii.gz: $(WORKDIR)/vesseltmp.nii.gz
	$(C3DEXE) -verbose $(WORKDIR)/anatomy.nii.gz  -thresh 140 inf 1 -1  $(WORKDIR)/hesstmp.nii.gz -replace 0 1 1 -1 -levelset-curvature 1.0 -levelset 10  -thresh -inf 0 1 0 -o $@

datalocation/radiomicsout.csv: datalocation/radiomics.csv
	python /opt/apps/pyradiomics/radiomics/scripts/commandlinebatch.py $< $@ -p Params.yaml  -l 1  -l 2 -v  5

