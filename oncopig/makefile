SHELL := /bin/bash
ANTSPATH=/opt/apps/ANTS/build/ANTS-build/Examples/
ATROPOSCMD=$(ANTSPATH)/Atropos -v 1 -d 3 -c [3,0.0]
# keep tmp files
.SECONDARY: 
vessels: otsu.1.nii.gz otsu.2.nii.gz otsu.3.nii.gz otsu.5.nii.gz

vessel.nii.gz:  otsu.1.nii.gz otsu.2.nii.gz otsu.3.nii.gz otsu.5.nii.gz 
	c3d -verbose $^ -accum -add -endaccum -binarize -dilate 1 3x3x1vox -erode 1 3x3x1vox -o $@
	c3d $@ -comp -dup -lstat > vessel.txt
	c3d $@ -comp -thresh 1 30 1 0  -o $@
otsu.%.nii.gz: vesselness.%.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
vesselness.5.nii.gz: sigmoid.nii.gz
	echo c3d -verbose $<  -hessobj 1 1.5 1.5 -o $@
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter $< $@ 1 1 .25 .25 10. 1.5
vesselness.%.nii.gz: sigmoid.nii.gz
	echo c3d -verbose $<  -hessobj 1 $* $* -o $@
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter $< $@ 1 1 .25 .25 10. $*

gmmPOSTERIORS4.nii.gz:   mask.nii.gz subtract.nii.gz speed.nii.gz
	$(ATROPOSCMD) -m [0.0,0x0x0] -i kmeans[4] -x $<  -a $(word 2,$^)  -a $(word 3,$^) -o [vessel.gmm.nii.gz,gmmPOSTERIORS%d.nii.gz] 
table.nii.gz: arterial.nii.gz 
	c3d -verbose $< -cmv -pop -thresh 100 410 1 0  -o $@
#	export HDF5_USE_FILE_LOCKING='FALSE'; python /rsrch1/ip/dtfuentes/github/livermask/applymodel.py --predictimage=/rsrch1/ip/dtfuentes/github/thermoembo/oncopig/arterial.nii.gz --segmentation=/rsrch1/ip/dtfuentes/github/thermoembo/oncopig/mask.nii.gz
mask.nii.gz: arterial.nii.gz table.nii.gz
	c3d -verbose  $< -thresh -900 inf 1 0 -erode 1 10x10x10vox -dilate 1 20x20x20vox $(word 2,$^) -multiply -type uchar -o $@
gradient.nii.gz: arterial.nii.gz 
	c3d -verbose $<  -smooth 1.2vox -grad -rms -o $@
sigmoid.nii.gz: subtract.nii.gz 
	c3d -verbose $<  -threshold -inf 95% 1 0 
	c3d -verbose $<  -scale -1 -shift 150.0 -scale 1 -exp -shift 1. -reciprocal -type uchar -o $@
subtract.nii.gz:  arterial.nii.gz noconresample.nii.gz 
	c3d -verbose  $^ -scale -1 -add -o $@


noconresample.nii.gz:  arterial.nii.gz nocontrast.nii.gz
	c3d $^ -reslice-identity -o $@
nocontrast.nii.gz:
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /mnt/FUS4/IPVL_research/ZPAF22S021/20220719/1.3.12.2.1107.5.1.4.83618.30000022071912215001100000001/1.3.12.2.1107.5.1.4.83618.30000022071912194835200000021/  $@ ;fi

arterial.nii.gz:
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /mnt/FUS4/IPVL_research/ZPAF22S021/20220719/1.3.12.2.1107.5.1.4.83618.30000022071912215001100000001/1.3.12.2.1107.5.1.4.83618.30000022071912194835200001346 $@ ;fi

