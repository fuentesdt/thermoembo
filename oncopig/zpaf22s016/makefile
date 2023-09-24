SHELL := /bin/bash
ANTSPATH=/opt/apps/ANTS/build/ANTS-build/Examples/
ATROPOSCMD=$(ANTSPATH)/Atropos -v 1 -d 3 -c 5
# keep tmp files
.SECONDARY: 

vessel.nii.gz:  otsu.1.nii.gz otsu.2.nii.gz otsu.3.nii.gz 
	c3d -verbose $^ -accum -add -endaccum -binarize -dilate 1 3x3x1vox -erode 1 3x3x1vox -o $@
	c3d -verbose manualregion.nii.gz $@ -overlap 1 > overlap.txt
	#c3d $@ -comp -dup -lstat > vessel.txt
	#c3d $@ -comp -thresh 1 30 1 0  -o $@
otsu.%.nii.gz: vesselness.%.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
vesselmax.nii.gz:  vesselness.1.nii.gz vesselness.2.nii.gz  vesselness.3.nii.gz
	c3d -verbose $^ -accum -max -endaccum -o $@
maxgmm: maxgmm2.nii.gz maxgmm3.nii.gz maxgmm4.nii.gz
maxgmm%.nii.gz: maskregion.nii.gz vesselmax.nii.gz 
	mkdir -p maxkmean$*
	$(ATROPOSCMD) -m [0.3,1x1x1] -i kmeans[$*] -x $<  -a $(word 2,$^)  -o [$@,maxkmean$*/gmmPOSTERIORS%d.nii.gz] 
otsumax.nii.gz: vesselmax.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
	c3d -verbose manualregion.nii.gz $@ -overlap 1 > overmax.txt
vesselnessns.%.nii.gz: 
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz         $@                1 1  .5    .5     5   $* 0
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz vesselnessnsgam.$*.nii.gz 1 1 1.e-9 1.e9    5   $* 0
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz vesselnessnsbet.$*.nii.gz 1 1 1.e-9  .5   1.e-9 $* 0
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz vesselnessnsalp.$*.nii.gz 1 1  .5   1.e9  1.e-9 $* 0
# note turn off the alpha component of the vesselness filter
vesselness.%.nii.gz: 
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz         $@              1 1 1.e-9 .5    5    $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz vesselnessgam.$*.nii.gz 1 1 1.e-9 1.e9  5    $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz vesselnessbet.$*.nii.gz 1 1 1.e-9 .5   1.e-9 $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz vesselnessalp.$*.nii.gz 1 1  .5   1.e9 1.e-9 $* 1
blobs: blobs.50.nii.gz blobs.100.nii.gz blobs.200.nii.gz
blobs.%.nii.gz: 
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz         $@         0 1 0.5   .5    5    $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz blobsgam.$*.nii.gz 0 1 1.e-9 1.e9  5    $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz blobsbet.$*.nii.gz 0 1 1.e-9 .5   1.e-9 $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz blobsalp.$*.nii.gz 0 1  .5   1.e9 1.e-9 $* 1
vesselnessnsgmm.%.nii.gz: 
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz            $@                1 1  .5   .5    5    $* 0
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz vesselnessnsgmmgam.$*.nii.gz 1 1 1.e-9 1.e9  5    $* 0
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz vesselnessnsgmmbet.$*.nii.gz 1 1 1.e-9 .5   1.e-9 $* 0
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz vesselnessnsgmmalp.$*.nii.gz 1 1   .5  1.e9 1.e-9 $* 0
vesselnessgmm.%.nii.gz: 
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz            $@              1 1  .5    .5   5    $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz vesselnessgmmgam.$*.nii.gz 1 1 1.e-9  1.e9 5    $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz vesselnessgmmbet.$*.nii.gz 1 1 1.e-9 .5   1.e-9 $* 1
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz vesselnessgmmalp.$*.nii.gz 1 1  .5   1.e9 1.e-9 $* 1
viewdbg:
	vglrun itksnap -g artregion.nii.gz     -o vesselness.1.nii.gz      vesselnessgam.1.nii.gz vesselnessbet.1.nii.gz vesselnessalp.1.nii.gz & vglrun itksnap -g artregion.nii.gz     -o vesselnessns.1.nii.gz  vesselnessnsgam.1.nii.gz  vesselnessnsbet.1.nii.gz  vesselnessnsalp.1.nii.gz & vglrun itksnap -g gmmpostregion.nii.gz -o vesselnessnsgmm.1.nii.gz vesselnessnsgmmgam.1.nii.gz vesselnessnsgmmbet.1.nii.gz vesselnessnsgmmalp.1.nii.gz & vglrun itksnap -g gmmpostregion.nii.gz -o vesselnessgmm.1.nii.gz   vesselnessgmmgam.1.nii.gz vesselnessgmmbet.1.nii.gz vesselnessgmmalp.1.nii.gz

cylinder1.nii.gz: manualregion.nii.gz
	c3d -verbose $<  -cmp  -pop -foreach -dup -times -endfor -accum -add -endaccum  -scale -.5 -exp  -scale  0.15915494309189535 -o $@
cylinderhess1scale.nii.gz: cylinder1.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter $< $@ 1 1 .5 .5 5 1 1
cylinderhess1noscale.nii.gz: cylinder1.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter $< $@ 1 1 .5 .5 5 1 0
manualregion.nii.gz: manual.nii.gz
	c3d $< -region 0x0x150vox 512x512x100vox -o $@
maskregion.nii.gz: mask.nii.gz
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
