


$(WORKDIR)/%/featureimages/HESSOBJ27.nii.gz: $(WORKDIR)/%/featureimages/maskfeatures
        for idvar in `seq 0 2`; do $(C3DEXE) -verbose $(WORKDIR)/$*/featureimages/Ven_DENOISE.nii.gz  -hessobj $$idvar .5 .5   -o $(WORKDIR)/$*/featureimages/HESSOBJ$${idvar}0.nii.gz; done
        for idvar in `seq 0 2`; do $(C3DEXE) -verbose $(WORKDIR)/$*/featureimages/Ven_DENOISE.nii.gz  -hessobj $$idvar  1  1   -o $(WORKDIR)/$*/featureimages/HESSOBJ$${idvar}1.nii.gz; done
        for idvar in `seq 0 2`; do $(C3DEXE) -verbose $(WORKDIR)/$*/featureimages/Ven_DENOISE.nii.gz  -hessobj $$idvar  2  2   -o $(WORKDIR)/$*/featureimages/HESSOBJ$${idvar}2.nii.gz; done
        for idvar in `seq 0 2`; do $(C3DEXE) -verbose $(WORKDIR)/$*/featureimages/Ven_DENOISE.nii.gz  -hessobj $$idvar  3  3   -o $(WORKDIR)/$*/featureimages/HESSOBJ$${idvar}3.nii.gz; done
        for idvar in `seq 0 2`; do $(C3DEXE) -verbose $(WORKDIR)/$*/featureimages/Ven_DENOISE.nii.gz  -hessobj $$idvar  4  4   -o $(WORKDIR)/$*/featureimages/HESSOBJ$${idvar}4.nii.gz; done
        for idvar in `seq 0 2`; do $(C3DEXE) -verbose $(WORKDIR)/$*/featureimages/Ven_DENOISE.nii.gz  -hessobj $$idvar  5  5   -o $(WORKDIR)/$*/featureimages/HESSOBJ$${idvar}5.nii.gz; done
        for idvar in `seq 0 2`; do $(C3DEXE) -verbose $(WORKDIR)/$*/featureimages/Ven_DENOISE.nii.gz  -hessobj $$idvar  6  6   -o $(WORKDIR)/$*/featureimages/HESSOBJ$${idvar}6.nii.gz; done
        for idvar in `seq 0 2`; do $(C3DEXE) -verbose $(WORKDIR)/$*/featureimages/Ven_DENOISE.nii.gz  -hessobj $$idvar  7  7   -o $(WORKDIR)/$*/featureimages/HESSOBJ$${idvar}7.nii.gz; done
$(WORKDIR)/%/vesselcomp.nii.gz: $(WORKDIR)/%/featureimages/HESSOBJ17.nii.gz
        $(C3DEXE) $< -thresh 1 inf 1 0 -comp  -o $@
