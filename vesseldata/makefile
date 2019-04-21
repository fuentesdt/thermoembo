ATROPOSCMD=$(ANTSPATH)/Atropos -d 3  -c [3,0.0] 
all: right.nii.gz left.nii.gz

# ex vivo
right.raw.nii.gz:
	DicomSeriesReadImageWrite2 /FUS4/IPVL_research/ZZPAF_LAYMAN_EXP42/20170420/1.3.12.2.1107.5.1.4.75774.30000017042115513085000030407/1.3.12.2.1107.5.1.4.75774.30000017042115513085000028413 $@
# ex vivo
left.raw.nii.gz:
	DicomSeriesReadImageWrite2 /FUS4/IPVL_research/ZZPAF_LAYMAN_EXP42/20170420/1.3.12.2.1107.5.1.4.75774.30000017042115513085000019365/1.3.12.2.1107.5.1.4.75774.30000017042115513085000018698 $@
right.roi.nii.gz:
	c3d right.raw.nii.gz -region 20x20x20vox 50x60x70vox -o $@
left.roi.nii.gz:
	c3d left.raw.nii.gz  -region 200x200x145vox 200x200x115vox -o $@

%.mask.nii.gz:  %.roi.nii.gz
	c3d -verbose $< -thresh -100 200 1 0 -comp -thresh 1 1 1 0 -dilate 1 3x3x3vox -erode 1 5x5x5vox -o $@

%.hess.nii.gz:  %.roi.nii.gz
	c3d -verbose $<  -hessobj 1  .5 2   -o  $@

%.gmm.nii.gz:   %.mask.nii.gz %.roi.nii.gz 
	$(ATROPOSCMD) -m [0.1,1x1x1] -i kmeans[3] -x $< -a $(word 2,$^) -o [$@,gmm$*POSTERIORS%d.nii.gz] 
	echo vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g $(word 2,$^) -s $@ -o  $(word 3,$^)

%.vessel.nii.gz: %.gmm.nii.gz
	c3d $< -thresh 1 1 1 0 $@
	vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g $*.roi.nii.gz -s $@ 
	echo vglrun /opt/apps/amira560/bin/start -tclcmd "load $*.roi.nii.gz; load $*.vessel.nii.gz;  create HxCastField; CastField data connect $*.vessel.nii.gz; CastField outputType setIndex 0 6; CastField create setLabel; Truth.nii.Labelfield ImageData connect $*.roi.nii.gz; "

landmarks:
	-c3d left.roi.nii.gz -scale 0 -type uchar left.lm.nii.gz
	-c3d ../temperaturedata/Kidney1Left_04202017_Exp42/temperature.0050.vtk -scale 0 -type uchar left.mrtilm.nii.gz
	vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g left.roi.nii.gz  -s left.lm.nii.gz & vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g ../temperaturedata/Kidney1Left_04202017_Exp42/temperature.0050.vtk -s left.mrtilm.nii.gz


%.rigid.mat:
	c3d -verbose $*.lm.nii.gz $*.mrtilm.nii.gz -alm 6 $@

left.mrti.nii.gz: left.rigid.mat
	c3d -verbose left.roi.nii.gz ../temperaturedata/Kidney1Left_04202017_Exp42/temperature.0050.vtk -reslice-matrix left.rigid.mat -o $@

viewreg:
	vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g left.roi.nii.gz  -o left.mrti.nii.gz
