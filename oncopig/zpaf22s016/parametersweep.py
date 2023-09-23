
import os

allradlist = [1,2,3,4,5,6,7,8,9]
radlistlist = [allradlist[:6], allradlist ]
gamlist = [5,50,500]
ablist  = [1,5,9]

for idab in ablist:
  for idgam in gamlist:
    for idrad in allradlist:
       nessfile = 'vesselness%d%d.%d.nii.gz' % (idab,idgam,idrad) 
       nesscmd = '/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter gmmpostregion.nii.gz %s  1 1 .%d .%d %f %d' % (nessfile,idab,idab,idgam*1.e-4,idrad)
       print(nesscmd)
       #if not os.path.isfile(nessfile):
       os.system(nesscmd)
       otsucmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter vesselness%d%d.%d.nii.gz otsu%d%d.%d.nii.gz 1  0' % (idab,idgam,idrad,idab,idgam,idrad)
       print(otsucmd)
       os.system(otsucmd)

for idab in ablist:
  for idgam in gamlist:
    for (idlist,radlist) in enumerate(radlistlist):
      nessniilist =  ' '.join(['vesselness%d%d.%d.nii.gz'%(idab,idgam,idrad) for idrad in radlist])
      maxcmd     = 'c3d -verbose %s  -accum -max -endaccum  -o vesselmax%d%d%d.nii.gz' %(nessniilist ,idab,idgam,idlist)
      print(maxcmd)
      os.system(maxcmd)
      otsumaxcmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter vesselmax%d%d%d.nii.gz otsumax%d%d%d.nii.gz 1  0' % (idab,idgam,idlist,idab,idgam,idlist)
      print(otsumaxcmd)
      os.system(otsumaxcmd)
      otsuniilist =  ' '.join(['otsu%d%d.%d.nii.gz'%(idab,idgam,idrad) for idrad in radlist])
      addbincmd = 'c3d -verbose %s -accum -add -endaccum -binarize -dilate 1 3x3x1vox -erode 1 3x3x1vox -o vessel%d%d%d.nii.gz' %(otsuniilist,idab,idgam,idlist)
      print(addbincmd)
      os.system(addbincmd)
      overbincmd ='c3d -verbose manualregion.nii.gz  vessel%d%d%d.nii.gz -overlap 1 > overlap%02d%03d%02d.txt' %(idab,idgam,idlist,idab,idgam,idlist)
      overmaxcmd ='c3d -verbose manualregion.nii.gz otsumax%d%d%d.nii.gz -overlap 1 > overmax%02d%03d%02d.txt' %(idab,idgam,idlist,idab,idgam,idlist)
      print(overbincmd )
      print(overmaxcmd )
      os.system(overbincmd )
      os.system(overmaxcmd )


