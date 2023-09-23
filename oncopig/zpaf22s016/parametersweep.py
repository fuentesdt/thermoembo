
import os

radlist = [1,2,3,4,5]
gamlist = [5,10,15,20]
ablist  = [1,5,9]

for idab in ablist:
  for idgam in gamlist:
    for idrad in radlist:
       nesscmd = '/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz vesselness%d%d.%d.nii.gz 1 1 .%d .%d %d %d' % (idab,idgam,idrad,idab,idab,idgam,idrad)
       print(nesscmd)
       os.system(nesscmd)
       otsucmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter vesselness%d%d.%d.nii.gz otsu%d%d.%d.nii.gz 1  0' % (idab,idgam,idrad,idab,idgam,idrad)
       print(otsucmd)
       os.system(otsucmd)

    maxcmd     = 'c3d -verbose vesselness%d%d.1.nii.gz vesselness%d%d.2.nii.gz vesselness%d%d.3.nii.gz vesselness%d%d.4.nii.gz vesselness%d%d.5.nii.gz -accum -max -endaccum -dilate 1 3x3x1vox -erode 1 3x3x1vox -o vesselmax%d%d.nii.gz' %(idab,idgam,idab,idgam,idab,idgam,idab,idgam,idab,idgam,idab,idgam)
    print(maxcmd)
    os.system(maxcmd)
    otsumaxcmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter vesselmax%d%d.nii.gz otsumax%d%d.nii.gz 1  0' % (idab,idgam,idab,idgam)
    print(otsumaxcmd)
    os.system(otsumaxcmd)
    addbincmd = 'c3d -verbose otsu%d%d.1.nii.gz otsu%d%d.2.nii.gz otsu%d%d.3.nii.gz otsu%d%d.4.nii.gz otsu%d%d.5.nii.gz -accum -add -endaccum -binarize -dilate 1 3x3x1vox -erode 1 3x3x1vox -o vessel%d%d.nii.gz' %(idab,idgam,idab,idgam,idab,idgam,idab,idgam,idab,idgam,idab,idgam)
    print(addbincmd)
    os.system(addbincmd)
    overbincmd ='c3d -verbose manualregion.nii.gz vessel%d%d.nii.gz  -overlap 1 > overlap%d%d.txt' %(idab,idgam,idab,idgam)
    overmaxcmd ='c3d -verbose manualregion.nii.gz otsumax%d%d.nii.gz -overlap 1 > overmax%d%d.txt' %(idab,idgam,idab,idgam)
    print(overbincmd )
    print(overmaxcmd )
    os.system(overbincmd )
    os.system(overmaxcmd )




