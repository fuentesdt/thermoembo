SHELL := /bin/bash
ANTSPATH=/opt/apps/ANTS/dev/install/bin/
ATROPOSCMD=$(ANTSPATH)/Atropos -v 1 -d 3 -c [3,0.0]
#SERIESLIST = $(shell sed 1d datalocation/wideformat.csv | cut -d, -f1 )
#READABLELIST = $(shell sed 1d datalocation/wideformat.csv | cut -d, -f2 )
#SERIESNUMLIST = $(shell sed 1d datalocation/wideformat.csv | cut -d, -f3 )
#PTUIDLIST = $(shell sed 1d datalocation/wideformat.csv | cut -d, -f4 )
#COUNT := $(words $(SERIESLIST))
#SEQUENCE = $(shell seq -f "%04g" $(COUNT))
#PROCESSLIST = 0003  0004   0005 0036
# keep tmp files
.SECONDARY: 
# ex vivo
vessel_centerlines.vtp: vesselmanual.vtk
	powershell.exe "C:\Users\dtfuentes\AppData\Local\Continuum\vmtk\python.exe C:\Users\dtfuentes\AppData\Local\Continuum\vmtk\Lib\site-packages\vmtk\vmtkcenterlines.py -ifile  $<  -ofile  $@"
# http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3
# c3d - command line image processing tools
# itksnap -  visualize  data
#

image.nii.gz:
	DicomSeriesReadImageWrite2 /mnt/FUS4/IPVL_research/ZPAF21S003/20210127/1.3.12.2.1107.5.1.4.83618.30000021012715452313500000001/1.3.12.2.1107.5.1.4.83618.30000021012714035248100000444/ $@
table.nii.gz: image.nii.gz
	c3d -verbose $< -cmv -pop -thresh -inf 440 1 0  -o table.nii.gz
	echo vglrun itksnap -g $< -s $@
mask.nii.gz:  image.nii.gz table.nii.gz
	c3d -verbose $< -thresh -900 inf 1 0 -erode 1 10x10x10vox -dilate 1 20x20x20vox   $(word 2,$^) -multiply -type uchar -o $@
	echo vglrun itksnap -g $< -s $@
gmm.nii.gz: image.nii.gz mask.nii.gz 
	$(ATROPOSCMD) -m [0,0x0x0] -i -i kmeans[3] -x $(word 2,$^)  -a $< -o [$@,$(@D)/gmmPOSTERIORS%d.nii.gz] 
	echo vglrun itksnap -g $< -s $@ -o gmmPOSTERIORS3.nii.gz
vesselness.%.nii.gz: gmmPOSTERIORS3.nii.gz 
	c3d -verbose $< -hessobj 1 $* $* -o $@
vesselness.5.nii.gz: gmmPOSTERIORS3.nii.gz 
	c3d -verbose $< -hessobj 1 .5 .5 -o $@
otsu.%.nii.gz: vesselness.%.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
vessel.nii.gz:  otsu.1.nii.gz otsu.2.nii.gz otsu.3.nii.gz otsu.4.nii.gz otsu.5.nii.gz 
	c3d -verbose $^ -accum -add -endaccum -binarize -comp -thresh 1 30 1 0 -dilate 1 3x3x1vox -erode 1 3x3x1vox -comp -o vessel.comp.nii.gz -thresh 1 1 1 0 -o $@
	echo vglrun itksnap -g image.nii.gz -s $@ 

Processed/i%/sigmoid.nii.gz: Processed/i%/image.nii.gz 
	c3d -verbose $< -dup -thresh -inf 500 1 0 -multiply -scale -1 -shift 120.0 -scale 1 -exp -shift 1. -reciprocal -type uchar -o $@

Processed/i%/mask.nii.gz: Processed/i%/otsusheet.5.nii.gz Processed/i%/otsusheet.6.nii.gz Processed/i%/otsusheet.7.nii.gz Processed/i%/bet_mask.nii.gz
	c3d -verbose  $< $(word 2,$^) -add $(word 3,$^) -add -binarize -replace 1 0 0 1 $(word 4,$^) -multiply -o $@

Processed/i%/sheet.1.nii.gz: Processed/i%/image.nii.gz 
	c3d -verbose $<  -hessobj 2 1 1  -o $@
Processed/i%/sheet.2.nii.gz: Processed/i%/image.nii.gz 
	c3d -verbose $<  -hessobj 2 2 2  -o $@
Processed/i%/sheet.3.nii.gz: Processed/i%/image.nii.gz 
	c3d -verbose $<  -hessobj 2 3 3  -o $@
Processed/i%/sheet.4.nii.gz: Processed/i%/image.nii.gz 
	c3d -verbose $<  -hessobj 2 4 4  -o $@
Processed/i%/sheet.5.nii.gz: Processed/i%/image.nii.gz 
	c3d -verbose $<  -hessobj 2 5 5  -o $@
Processed/i%/sheet.6.nii.gz: Processed/i%/image.nii.gz 
	c3d -verbose $<  -hessobj 2 6 6  -o $@
Processed/i%/sheet.7.nii.gz: Processed/i%/image.nii.gz 
	c3d -verbose $<  -hessobj 2 7 7  -o $@
Processed/i%/otsusheet.1.nii.gz: Processed/i%/sheet.1.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
Processed/i%/otsusheet.2.nii.gz: Processed/i%/sheet.2.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
Processed/i%/otsusheet.3.nii.gz: Processed/i%/sheet.3.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
Processed/i%/otsusheet.4.nii.gz: Processed/i%/sheet.4.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
Processed/i%/otsusheet.5.nii.gz: Processed/i%/sheet.5.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
Processed/i%/otsusheet.6.nii.gz: Processed/i%/sheet.6.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 
Processed/i%/otsusheet.7.nii.gz: Processed/i%/sheet.7.nii.gz 
	/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter $< $@ 1  0 

rfmask: $(addprefix Processed/i,$(addsuffix /rfmask.nii.gz,$(PROCESSLIST))) $(addprefix Processed/i,$(addsuffix /rfsamples.nii.gz,$(PROCESSLIST))) 
Processed/i%/rfmask.nii.gz: Processed/i%/image.nii.gz Processed/i%/bet_mask.nii.gz
	c3d -verbose $^  -erode 1 4x4x0vox -multiply  -o $@
Processed/i%/rfsamples.nii.gz: Processed/i%/gmm.nii.gz Processed/i%/bet_mask.nii.gz
	c3d -verbose $^  -erode 1 4x4x0vox -multiply  -o $@

view: $(addprefix Processed/i,$(addsuffix /view,$(PROCESSLIST)))
Processed/i%/view:
	cat $(@D)/header.json
	vglrun itksnap -g $(@D)/image.nii.gz -s $(@D)/gmm.nii.gz -o $(@D)/vesselness.?.nii.gz $(@D)/otsu.?.nii.gz $(@D)/otsusheet.?.nii.gz & vglrun itksnap -g $(@D)/rfmask.nii.gz  -s $(@D)/rfsamples.nii.gz 

subtract.nii.gz:
	c3d -verbose  Brats18_2013_13_1_t1ce.nii.gz Brats18_2013_13_1_t1.nii.gz -scale -1 -add -o $@
arclengthfiducials.nii.gz: subtract.nii.gz
	if [ ! -f $@  ] ; then c3d $< -scale 0 -type uchar $@ ; else touch $@ ; fi
arclength.json: sigmoid.nii.gz arclengthfiducials.nii.gz
	python PathExtraction.py $(basename $@)  $^ 

vessel.hess.nii.gz: 
	c3d -verbose  vesselness.1.nii.gz  -thresh 30 inf 1 0 -o $@

thin.nii.gz: vessel.nii.gz
	/rsrch3/ip/dtfuentes/github/ctperfusion/ThinImage $< $@ 
centerline.nii.gz: vessel.nii.gz thin.nii.gz 
	c3d $^ -binarize -add -o  $@
sagittal.nii.gz: 
	/rsrch3/ip/dtfuentes/github/ctperfusion/ConnectedThresholdImageFilter subtract.nii.gz  $@ 112 203 78 300 1000
# given a 3d vessel segmentation, create skeleton of vessel:
# 
#      https://itk.org/Doxygen/html/classitk_1_1BinaryThinningImageFilter.html
