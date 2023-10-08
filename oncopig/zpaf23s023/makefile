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
maxgmm: maxgmm2.nii.gz maxgmm3.nii.gz maxgmm4.nii.gz
maxgmm%.nii.gz: maskregion.nii.gz vesselmax.nii.gz 
	mkdir -p maxkmean$*
	$(ATROPOSCMD) -m [0.3,1x1x1] -i kmeans[$*] -x $<  -a $(word 2,$^)  -o [$@,maxkmean$*/gmmPOSTERIORS%d.nii.gz] 
otsuliver.nii.gz: vesselmax.nii.gz liver.nii.gz
	c3d -verbose $^ -times -o $@
otsucomp.nii.gz: 
	c3d -verbose otsumaxall0.nii.gz -as SEG -cmv -pop -pop  -thresh 0 205 1 0 -as MASK -push SEG -times -comp -thresh 1 5 1 0 -comp -as LEFT -push MASK -replace 1 0 0 1  -push SEG -times -comp -thresh 1 5 1 0 -comp -push LEFT -add -o $@

otsumax.nii.gz: otsuliver.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
	c3d $@ -dilate 1 3x3x1vox -erode 1 3x3x1vox -verbose -o $@
	c3d -verbose manual.nii.gz $@ -overlap 1 > overmax.txt
hessmag: hessianmagnitude.1.nii.gz hessianmagnitude.2.nii.gz hessianmagnitude.3.nii.gz
hessianmagnitude.%.nii.gz: artregion.nii.gz
	c3d -verbose $< -hesseig $* -oo eig$*%02d.nii.gz -foreach -dup -times -endfor -accum -add -endaccum -sqrt -o $@ 
manualrf.nii.gz: manualregion.nii.gz
	c3d -verbose $< -shift 1 -type uchar -o $@
myforest1.rf: vesselnessalp.1.nii.gz vesselnessalp.2.nii.gz vesselnessalp.3.nii.gz vesselnessbet.1.nii.gz vesselnessbet.2.nii.gz vesselnessbet.3.nii.gz vesselnessgam.1.nii.gz vesselnessgam.2.nii.gz vesselnessgam.3.nii.gz  manualrf.nii.gz
	c3d -verbose $^  -rf-train $@
class_prob.nii.gz: myforest1.rf
	c3d -verbose vesselnessalp.1.nii.gz vesselnessalp.2.nii.gz vesselnessalp.3.nii.gz vesselnessbet.1.nii.gz vesselnessbet.2.nii.gz vesselnessbet.3.nii.gz vesselnessgam.1.nii.gz vesselnessgam.2.nii.gz vesselnessgam.3.nii.gz  -rf-apply $<  -omc class_prob.nii.gz

rfotsu.nii.gz: class_prob.nii.gz
	c3d -verbose -mcs $< -oo class_prob%02d.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter class_prob01.nii.gz  $@ 1  0 
	c3d $@ -dilate 1 3x3x1vox -erode 1 3x3x1vox -verbose -o $@
	c3d -verbose manualregion.nii.gz $@ -overlap 1 > overrf.txt
# hessian magnitude inside liver should be used for vesselness gamma parameter 
gamma.%: hessianmagnitude.%.nii.gz
	c3d $< liverregion.nii.gz -lstat
viewdbg:
	vglrun itksnap -g artregion.nii.gz     -o vesselness.1.nii.gz      vesselnessgam.1.nii.gz vesselnessbet.1.nii.gz vesselnessalp.1.nii.gz & vglrun itksnap -g artregion.nii.gz     -o vesselnessns.1.nii.gz  vesselnessnsgam.1.nii.gz  vesselnessnsbet.1.nii.gz  vesselnessnsalp.1.nii.gz & vglrun itksnap -g gmmpostregion.nii.gz -o vesselnessnsgmm.1.nii.gz vesselnessnsgmmgam.1.nii.gz vesselnessnsgmmbet.1.nii.gz vesselnessnsgmmalp.1.nii.gz & vglrun itksnap -g gmmpostregion.nii.gz -o vesselnessgmm.1.nii.gz   vesselnessgmmgam.1.nii.gz vesselnessgmmbet.1.nii.gz vesselnessgmmalp.1.nii.gz

regressiontest: regressionhessmagnitude.nii.gz regressionotsu.nii.gz regressionvesselnessalpha.nii.gz regressionvesselnessbeta.nii.gz regressionvesselnessgamma.nii.gz
	vglrun itksnap -g regressionscale.nii.gz -s regressiontest.nii.gz -o regressioneig0?.nii.gz regressionvesselness*.nii.gz regressionhessmagnitude.nii.gz   & vglrun itksnap -g regressionscale.nii.gz -s regressionotsu.nii.gz -o regressioneig0?.nii.gz regressionvesselness*.nii.gz regressionhessmagnitude.nii.gz  
regressionscale.nii.gz: regressiontest.nii.gz
	c3d $< -binarize -scale 500 -type float -o $@
regressionhessmagnitude.nii.gz: regressionscale.nii.gz
	c3d -verbose $< -hesseig .818 -oo regressioneig%02d.nii.gz -foreach -dup -times -endfor -accum -add -endaccum -sqrt -o $@ 
regressionvesselnessalpha.nii.gz: regressionscale.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter $< $@ 1 1  .5   1.e9   0.  .818 1
regressionvesselnessbeta.nii.gz:  regressionscale.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter $< $@ 1 1 0.  .5    0.  .818 1
regressionvesselnessgamma.nii.gz: regressionscale.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter $< $@ 1 1 0. 1.e9 167.     .818 1
regressionvesselness.nii.gz: regressionscale.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter $< $@ 1 1  .5    .5  167.     .818 1
regressionotsu.nii.gz: regressionvesselness.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 

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
artliver.nii.gz: artregion.nii.gz liverregion.nii.gz
	c3d $^ -times -o $@
liverregion.nii.gz: liver.nii.gz
	c3d $< -region 0x0x150vox 512x512x100vox -o $@
gmmpostregion.nii.gz: kmean3/gmmPOSTERIORS3.nii.gz
	c3d $< -region 0x0x150vox 512x512x100vox -o $@

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

vessel: art.vessel.nii.gz ven.vessel.nii.gz
%.vessel.nii.gz: %.vesselmax.nii.gz 
	c3d -verbose $< -threshold 1.0 inf 1 0 -o $@
%.vesselmax.nii.gz:  %.vesselness.1.nii.gz %.vesselness.2.nii.gz  %.vesselness.3.nii.gz 
	c3d -verbose $^ -accum -max -endaccum -o $@
%.denoise.nii.gz: %.reg.nii.gz
	docker run --entrypoint=/opt/mitk/ImageDenoising/TotalVariationDenoisingImage  --rm -it --user $$(id -u):$$(id -g) -v $(PWD)/:/data/  -v $(PWD):/out iptools:latest /data/$< /out/$@ .1 30
pixelsize=0.90625
art.vesselness.%.nii.gz:
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter art.denoise.nii.gz $@ 1 1  .5   5.    5    $(shell python -c "print $* * $(pixelsize)" ) 1
art.vesselness.5.nii.gz:                                                                                                     
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter art.denoise.nii.gz $@ 1 1  .5   5.    5    $(shell python -c "print .5 * $(pixelsize)" ) 1
ven.vesselness.%.nii.gz:                                                                                                     
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter ven.denoise.nii.gz $@ 1 1  .5   5.    5    $(shell python -c "print $* * $(pixelsize)" ) 1
ven.vesselness.5.nii.gz:                                                                                                     
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter ven.denoise.nii.gz $@ 1 1  .5   5.    5    $(shell python -c "print .5 * $(pixelsize)" ) 1
subtract:  art.subtract.nii.gz ven.subtract.nii.gz
%.subtract.nii.gz:  %.reg.nii.gz pre.reg.nii.gz 
	c3d -verbose  $^ -scale -1 -add -o $@

art.lung.nii.gz: art.body.nii.gz
	c3d arterial.nii.gz artbody.nii.gz -times -o bodyroi.nii.gz
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter bodyroi.nii.gz  $@ 1  0 
	c3d $@ -replace 1 0 0 1 -comp -o $@ 

reg: pre.reg.nii.gz ven.reg.nii.gz  art.reg.nii.gz
art.reg.nii.gz: art.nii.gz
	cp $< $@
%.reg.nii.gz: %.nii.gz  %.body.nii.gz art.nii.gz  art.body.nii.gz 
	export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=28; /opt/apps/ANTS/build/ANTS-build/Examples/antsRegistration --verbose 1 --dimensionality 3 --float 0 --collapse-output-transforms 1 --output [$(basename $(basename $@)),$@] --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] -x [art.body.nii.gz,$*.body.nii.gz] -r [ art.body.nii.gz,$*.body.nii.gz,1] --transform Rigid[ 0.1 ] --metric MI[ art.nii.gz,$<,1,32,Regular,0.25 ] --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[ 0.1 ] --metric MI[ art.nii.gz,$<,1,32,Regular,0.25 ] --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[ 0.1,3,0 ] --metric CC[art.nii.gz,$<,1,4 ] --convergence [ 100x70x50x20,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox 
body: pre.body.nii.gz art.body.nii.gz ven.body.nii.gz
%.body.nii.gz: %.nii.gz
	c3d -verbose $< -as REG -threshold -475 inf 1 0 -type char -pad 0x0x1vox 0x0x1vox 1  -holefill 1 0 -popas MASK -push REG -push MASK -reslice-identity -o $*body.1.nii.gz
	c3d -verbose $< -as REG -threshold -190 inf 1 0 -type char -pad 0x0x1vox 0x0x1vox 1  -holefill 1 0 -popas MASK -push REG -push MASK -reslice-identity -o $*body.2.nii.gz
	c3d -verbose $< -as REG -threshold -165 inf 1 0 -type char -pad 0x0x1vox 0x0x1vox 1  -holefill 1 0 -popas MASK -push REG -push MASK -reslice-identity -o $*body.3.nii.gz
	c3d -verbose $< -as REG -threshold -175 inf 1 0 -type char -pad 0x0x1vox 0x0x1vox 1  -holefill 1 0 -popas MASK -push REG -push MASK -reslice-identity -o $*body.4.nii.gz
	c3d -verbose $< -as REG -threshold -150 inf 1 0 -type char -pad 0x0x1vox 0x0x1vox 1  -holefill 1 0 -popas MASK -push REG -push MASK -reslice-identity -o $*body.5.nii.gz
	c3d -verbose $< -as REG -threshold -155 inf 1 0 -type char -pad 0x0x1vox 0x0x1vox 1  -holefill 1 0 -popas MASK -push REG -push MASK -reslice-identity -o $*body.6.nii.gz
	c3d -verbose $< -as REG -threshold -400 inf 1 0 -type char -pad 0x0x1vox 0x0x1vox 1  -holefill 1 0 -popas MASK -push REG -push MASK -reslice-identity -o $*body.7.nii.gz
	c3d -verbose $< -as REG -threshold -200 inf 1 0 -type char -pad 0x0x1vox 0x0x1vox 1  -holefill 1 0 -popas MASK -push REG -push MASK -reslice-identity -o $*body.8.nii.gz
	c3d -verbose $*body.?.nii.gz -accum -add -endaccum -thresh 8 8 1 0 -comp -thresh 1 1 1 0 -type char -o $@

raw: nocontrast.nii.gz arterial.nii.gz venous.nii.gz 
resample: art.nii.gz pre.resample.nii.gz ven.resample.nii.gz 
	vglrun itksnap -g $> -o $(word 2,$^) $(word 3,$^)
info:
	c3d nocontrast.nii.gz -info  arterial.nii.gz  -info venous.nii.gz  -info lateart.nii.gz    -info 
%.resample.nii.gz:  art.nii.gz %.nii.gz
	c3d $^ -reslice-identity -o $@
pre.nii.gz:
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /rsrch1/ip/ip-comp_rsch_lab/oncopigdata/ZPAF23S023/20230510/1.3.12.2.1107.5.1.4.83618.30000023051011532834500000001/1.3.12.2.1107.5.1.4.83618.30000023050916382426000007919 $@ ;fi

art.nii.gz:
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /rsrch1/ip/ip-comp_rsch_lab/oncopigdata/ZPAF23S023/20230510/1.3.12.2.1107.5.1.4.83618.30000023051011532834500000001/1.3.12.2.1107.5.1.4.83618.30000023050916382426000013350 $@ ;fi

ven.nii.gz:
	if [ ! -f $@  ] ; then DicomSeriesReadImageWrite2 /rsrch1/ip/ip-comp_rsch_lab/oncopigdata/ZPAF23S023/20230510/1.3.12.2.1107.5.1.4.83618.30000023051011532834500000001/1.3.12.2.1107.5.1.4.83618.30000023050916382426000011535 $@ ;fi
