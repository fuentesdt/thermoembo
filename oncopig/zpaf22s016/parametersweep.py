
import os

radlist = [1,2,3,4,5,6,7,8,9]
gamlist = [5,10,15,20]
ablist  = [1,5,9]

for idab in ablist:
  for idgam in gamlist:
    for idrad in radlist:
       nessfile = 'vesselness%d%d.%d.nii.gz' % (idab,idgam,idrad) 
       nesscmd = '/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz %s  1 1 .%d .%d %d %d' % (nessfile,idab,idab,idgam,idrad)
       print(nesscmd)
       if not os.path.isfile(nessfile):
         os.system(nesscmd)
       otsucmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter vesselness%d%d.%d.nii.gz otsu%d%d.%d.nii.gz 1  0' % (idab,idgam,idrad,idab,idgam,idrad)
       print(otsucmd)
       os.system(otsucmd)

    nessniilist =  ' '.join(['vesselness%d%d.%d.nii.gz'%(idab,idgam,idrad) for idrad in radlist])
    maxcmd     = 'c3d -verbose %s  -accum -max -endaccum  -o vesselmax%d%d.nii.gz' %(nessniilist ,idab,idgam)
    print(maxcmd)
    os.system(maxcmd)
    otsumaxcmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter vesselmax%d%d.nii.gz otsumax%d%d.nii.gz 1  0' % (idab,idgam,idab,idgam)
    print(otsumaxcmd)
    os.system(otsumaxcmd)
    otsuniilist =  ' '.join(['otsu%d%d.%d.nii.gz'%(idab,idgam,idrad) for idrad in radlist])
    addbincmd = 'c3d -verbose %s -accum -add -endaccum -binarize -dilate 1 3x3x1vox -erode 1 3x3x1vox -o vessel%d%d.nii.gz' %(otsuniilist,idab,idgam)
    print(addbincmd)
    os.system(addbincmd)
    overbincmd ='c3d -verbose manualregion.nii.gz vessel%d%d.nii.gz  -overlap 1 > overlap%02d%02d.txt' %(idab,idgam,idab,idgam)
    overmaxcmd ='c3d -verbose manualregion.nii.gz otsumax%d%d.nii.gz -overlap 1 > overmax%02d%02d.txt' %(idab,idgam,idab,idgam)
    print(overbincmd )
    print(overmaxcmd )
    os.system(overbincmd )
    os.system(overmaxcmd )




