SHELL := /bin/bash
ANTSPATH=/opt/apps/ANTS/build/ANTS-build/Examples/
ATROPOSCMD=$(ANTSPATH)/Atropos -v 1 -d 3 -c 5
Gamma=10
# keep tmp files
.SECONDARY: 
vessels: otsu.1.nii.gz otsu.2.nii.gz otsu.3.nii.gz otsu.4.nii.gz otsu.5.nii.gz

vessel.nii.gz:  otsu.1.nii.gz otsu.2.nii.gz otsu.3.nii.gz otsu.5.nii.gz 
	c3d -verbose $^ -accum -add -endaccum -binarize -dilate 1 3x3x1vox -erode 1 3x3x1vox -o $@
	c3d -verbose manualregion.nii.gz $@ -overlap 1 > overlap$(Gamma).txt
	#c3d $@ -comp -dup -lstat > vessel.txt
	#c3d $@ -comp -thresh 1 30 1 0  -o $@
otsu.%.nii.gz: vesselness.%.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
vesselmax.nii.gz:  vesselness.1.nii.gz vesselness.2.nii.gz vesselness.3.nii.gz vesselness.5.nii.gz 
	c3d -verbose $^ -accum -max -endaccum -dilate 1 3x3x1vox -erode 1 3x3x1vox -o $@
otsumax.nii.gz: vesselmax.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
	c3d -verbose manualregion.nii.gz $@ -overlap 1 > overmax$(Gamma).txt
vesselness.%.nii.gz: 
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz $@ 1 1 .5 .5 $(Gamma) $*

manualregion.nii.gz: manual.nii.gz
	c3d $< -region 0x0x150vox 512x512x100vox -o $@
artregion.nii.gz: arterial.nii.gz
	c3d $< -region 0x0x150vox 512x512x100vox -o $@
gmmpostregion.nii.gz: kmean3/gmmPOSTERIORS3.nii.gz
	c3d $< -region 0x0x150vox 512x512x100vox -o $@
hessianmagnitude.%.nii.gz: gmmpostregion.nii.gz
	c3d -verbose $< -hesseig $* -oo eig$*%02d.nii.gz -foreach -dup -times -endfor -accum -add -endaccum -sqrt -o $@ 

gmmthresh.nii.gz: vessel.gmm.nii.gz
	c3d $< -thresh 3 4 1 0  -o $@
gmmoverlap:
	c3d vesselgmm2.nii.gz manual.nii.gz              -overlap 1
	c3d vesselgmm2.nii.gz manual.nii.gz -replace 1 2 -overlap 2
	c3d vesselgmm3.nii.gz manual.nii.gz              -overlap 1
	c3d vesselgmm3.nii.gz manual.nii.gz -replace 1 2 -overlap 2
	c3d vesselgmm3.nii.gz manual.nii.gz -replace 1 3 -overlap 3
	c3d vesselgmm4.nii.gz manual.nii.gz              -overlap 1
	c3d vesselgmm4.nii.gz manual.nii.gz -replace 1 2 -overlap 2
	c3d vesselgmm4.nii.gz manual.nii.gz -replace 1 3 -overlap 3
	c3d vesselgmm4.nii.gz manual.nii.gz -replace 1 4 -overlap 4
	c3d vesselgmm5.nii.gz manual.nii.gz              -overlap 1
	c3d vesselgmm5.nii.gz manual.nii.gz -replace 1 2 -overlap 2
	c3d vesselgmm5.nii.gz manual.nii.gz -replace 1 3 -overlap 3
	c3d vesselgmm5.nii.gz manual.nii.gz -replace 1 4 -overlap 4
	c3d vesselgmm5.nii.gz manual.nii.gz -replace 1 5 -overlap 5
vesselgmm: vesselgmm2.nii.gz vesselgmm3.nii.gz vesselgmm4.nii.gz vesselgmm5.nii.gz
vesselgmm%.nii.gz: mask.nii.gz arterial.nii.gz 
	mkdir -p kmean$*
	$(ATROPOSCMD) -m [0.3,1x1x1] -i kmeans[$*] -x $<  -a $(word 2,$^)  -o [vesselgmm$*.nii.gz,kmean$*/gmmPOSTERIORS%d.nii.gz] 
table.nii.gz: arterial.nii.gz 
	c3d -verbose $< -cmv  -thresh  150 230  1 0 -popas ZDIR -thresh 110 480 1 0 -push ZDIR -times -o $@
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


raw: nocontrast.nii.gz arterial.nii.gz venous.nii.gz lateart.nii.gz
resample: noconresample.nii.gz venresample.nii.gz lateresample.nii.gz
info:
	c3d nocontrast.nii.gz -info  arterial.nii.gz  -info venous.nii.gz  -info lateart.nii.gz    -info 
noconresample.nii.gz:  arterial.nii.gz nocontrast.nii.gz
	c3d $^ -reslice-identity -o $@
venresample.nii.gz:  arterial.nii.gz venous.nii.gz
	c3d $^ -reslice-identity -o $@
lateresample.nii.gz:  arterial.nii.gz lateart.nii.gz
	c3d $^ -reslice-identity -o $@
nocontrast.nii.gz:
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /rsrch1/ip/ip-comp_rsch_lab/oncopigdata/ZPAF22S016/20220809/1.3.12.2.1107.5.1.4.83618.30000022080912274636100000001/1.3.12.2.1107.5.1.4.83618.30000022080912225286100000021 $@ ;fi
arterial.nii.gz:
	#if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /mnt/FUS4/IPVL_research/ZPAF22S016/20220809/1.3.12.2.1107.5.1.4.83618.30000022080912274636100000001/1.3.12.2.1107.5.1.4.83618.30000022080912225286100001398/  $@ ;fi
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /rsrch1/ip/ip-comp_rsch_lab/oncopigdata/ZPAF22S016/20220809/1.3.12.2.1107.5.1.4.83618.30000022080912274636100000001/1.3.12.2.1107.5.1.4.83618.30000022080912225286100001398 $@ ;fi
venous.nii.gz:
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /rsrch1/ip/ip-comp_rsch_lab/oncopigdata/ZPAF22S016/20220809/1.3.12.2.1107.5.1.4.83618.30000022080912274636100000001/1.3.12.2.1107.5.1.4.83618.30000022080912225286100000939 $@ ;fi
lateart.nii.gz:
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /rsrch1/ip/ip-comp_rsch_lab/oncopigdata/ZPAF22S016/20220809/1.3.12.2.1107.5.1.4.83618.30000022080912274636100000001/1.3.12.2.1107.5.1.4.83618.30000022080912225286100001856 $@ ;fi
