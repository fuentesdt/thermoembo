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
maxgmm%.nii.gz: liverregion.nii.gz artregion.nii.gz 
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
gmmpre.nii.gz: maxgmm2.nii.gz
	c3d $< -thresh 2 2 1 0 -o $@
gmmvesselness.%.nii.gz: 
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpre.nii.gz $@ 1 1  .5   .5    5    $* 1
gmmvessel.nii.gz:  gmmvesselness.1.nii.gz gmmvesselness.2.nii.gz  gmmvesselness.3.nii.gz gmmvesselness.4.nii.gz  gmmvesselness.5.nii.gz   gmmvesselness.6.nii.gz  gmmvesselness.7.nii.gz   gmmvesselness.8.nii.gz  gmmvesselness.9.nii.gz gmmvesselness.10.nii.gz
	c3d -verbose $^ -accum -max -endaccum -o $@
gmmotsumax.nii.gz: gmmvessel.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
gmmthin.nii.gz: gmmotsumax.nii.gz
	/rsrch3/ip/dtfuentes/github/ctperfusion/ThinImage $< $@
soi.nii.gz: gmmthin.nii.gz
	c3d $< -dilate 255 4x4x4vox  -o $@
centerline.vtk:
	export PYTHONPATH=/rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/lib/python2.7/site-packages;
	/rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtkmarchingcubes -l 1 -ifile gmmotsumax.nii.gz -ofile model.vtp 
	vglrun /rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtkcenterlines -ifile model.vtp -seedselector 'idlist' -targetids 17466 -sourceids 16143 -ofile foo_centerlines.vtp
	vglrun /rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtklevelsetsegmentation -ifile gmmotsumax.nii.gz --pipe /rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtkmarchingcubes -i @.o -ofile model.vtp
	/rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtkmarchingcubes -l 1 -ifile regressionvmtk.nii.gz -ofile regmodel.vtp 
	/rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtksurfacedecimation -ifile regmodel.vtp  -ofile decmodel.vtp  -reduction .9
	/rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtksurfacesmoothing -ifile decmodel.vtp  -ofile smtmodel.vtp  -iterations 10
	vglrun /rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtkcenterlines -ifile regmodel.vtp -ofile reg_centerlines.vtp
	vglrun /rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtkcenterlines -ifile gmmotsumax.nii.surf.stl -ofile acenterlines.vtp
centerlines.stl:
	export AMIRA_DATADIR=$(PWD);vglrun /opt/apps/Amira/2020.2/bin/start -tclcmd 'load gmmotsumax.nii.gz; create HxCastField ConvertImage; ConvertImage data connect gmmotsumax.nii.gz; ConvertImage fire; ConvertImage outputType setIndex 0 7; ConvertImage create result setLabel ; gmmotsumax.nii.to-labelfield-8_bits ImageData connect gmmotsumax.nii.gz; set hideNewModules 0 ; create HxGMC GenerateSurface ; GenerateSurface setIconPosition 399 68 ; GenerateSurface setVar CustomHelp {HxGMC} ; GenerateSurface setVar CustomXML {HxGMC_config.xml} ; GenerateSurface data connect gmmotsumax.nii.to-labelfield-8_bits ; GenerateSurface fire ; GenerateSurface smoothing setIndex 0 3 ; GenerateSurface smoothingExtent setMinMax 1 9 ; GenerateSurface smoothingExtent setButtons 0 ; GenerateSurface smoothingExtent setEditButton 1 ; GenerateSurface smoothingExtent setIncrement 0.533333 ; GenerateSurface smoothingExtent setValue 2 ; GenerateSurface smoothingExtent setSubMinMax 1 9 ; GenerateSurface options setState {item 0 0 item 1 0 } ; GenerateSurface borderOnOff setValue 1 ; GenerateSurface algorithmMode setIndex 0 0 ; GenerateSurface borderSettings setValue 0 1 ; GenerateSurface borderSettings setToggleVisible 0 1 ; GenerateSurface borderSettings setValue 1 0 ; GenerateSurface borderSettings setToggleVisible 1 1 ; GenerateSurface borderSettings setValue 2 0 ; GenerateSurface borderSettings setToggleVisible 2 1 ; GenerateSurface borderSettings setValue 3 0 ; GenerateSurface borderSettings setToggleVisible 3 0 ; GenerateSurface materialList setIndex 0 0 ; GenerateSurface smoothMaterial setIndex 0 0 ; GenerateSurface applyTransformToResult 1 ; GenerateSurface fire ; GenerateSurface setViewerMask 16383 ; GenerateSurface select ; GenerateSurface setPickable 1 ; set hideNewModules 0 ; [ GenerateSurface create gmmotsumax.nii.surf ] setLabel gmmotsumax.nii.surf ; gmmotsumax.nii.surf setIconPosition 20 115 ; gmmotsumax.nii.surf master connect GenerateSurface result 0 ; gmmotsumax.nii.surf setTransform -1 0 0 0 0 -1 0 0 0 0 1 0 225.091 304.091 464.5 1 ; gmmotsumax.nii.surf fire ; gmmotsumax.nii.surf setViewerMask 16383; gmmotsumax.nii.surf exportData "STL ascii" $@';

centerlines.vtp: centerlines.stl
	export PYTHONPATH=/rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/lib/python2.7/site-packages; vglrun python /rsrch3/ip/dtfuentes/github/vmtk/vmtk-build/Install/bin/vmtkcenterlines -ifile $< -ofile $@


# hessian magnitude inside liver should be used for vesselness gamma parameter 
gamma.%: hessianmagnitude.%.nii.gz
	c3d $< liverregion.nii.gz -lstat
vesselness.%.nii.gz: 
	/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter arterial.nii.gz $@ 1 1  .5   .5    5    $* 1
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

dsa.nii.gz:
	DicomSeriesReadImageWrite2 /mnt/FUS4/IPVL_research/ZPAF22S016/20220809/1.3.12.2.1107.5.4.5.109411.30000022080912120636200000003/1.3.12.2.1107.5.4.5.109411.30000022080912315231600000033/ $@

dsaroi.nii.gz: dsa.nii.gz
	c3d $< -cmv -thresh 6 25 1 0 $< -thresh 10 inf 1 0 -multiply -o $@
dsagmm.nii.gz: dsaroi.nii.gz dsa.nii.gz 
	$(ATROPOSCMD) -m [0.3,1x1x1] -i kmeans[2] -x $<  -a $(word 2,$^)  -o [$@,gmmPOSTERIORS%d.nii.gz] 
