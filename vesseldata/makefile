ATROPOSCMD=$(ANTSPATH)/Atropos -d 3  -c [3,0.0] 
all: right.nii.gz left.nii.gz

# ex vivo
right.nii.gz:
	DicomSeriesReadImageWrite2 /FUS4/IPVL_research/ZZPAF_LAYMAN_EXP42/20170420/1.3.12.2.1107.5.1.4.75774.30000017042115513085000030407/1.3.12.2.1107.5.1.4.75774.30000017042115513085000028413 $@
# ex vivo
left.nii.gz:
	DicomSeriesReadImageWrite2 /FUS4/IPVL_research/ZZPAF_LAYMAN_EXP42/20170420/1.3.12.2.1107.5.1.4.75774.30000017042115513085000019365/1.3.12.2.1107.5.1.4.75774.30000017042115513085000018698 $@
maskleft.nii.gz: 
	c3d -verbose left.nii.gz -thresh -100 200 1 0 -comp -thresh 1 1 1 0 -dilate 1 3x3x3vox -erode 1 3x3x3vox -o $@
maskright.nii.gz: 
	c3d -verbose right.nii.gz -thresh -100 200 1 0 -comp -thresh 1 1 1 0 -dilate 1 3x3x3vox -erode 1 3x3x3vox -o $@
gmmleft.nii.gz: 
	$(ATROPOSCMD) -m [0.1,1x1x1] -i kmeans[3] -x maskleft.nii.gz -a left.nii.gz  -o [$@,gmmleftPOSTERIORS%d.nii.gz] 
	echo vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g left.nii.gz -s gmmleft.nii.gz
gmmright.nii.gz: 
	$(ATROPOSCMD) -m [0.1,1x1x1] -i kmeans[3] -x maskright.nii.gz -a right.nii.gz  -o [$@,gmmrightPOSTERIORS%d.nii.gz] 
	echo vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap -g right.nii.gz -s gmmright.nii.gz
