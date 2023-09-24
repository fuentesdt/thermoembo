
import os

allradlist = [.4,.8,1.,1.6,2.0,2.4,3.]
radlistlist = [allradlist[2::2], allradlist ]
gamlist = [5,50,500]
alplist  = [1,5,9]
betlist  = [1,5,9]
paramlist ={'alp':(.5   ,1.e9,1.e-9),
            'bet':(1.e-9, .5 ,1.e-9),
            'gam':(1.e-9,1.e9,5.   ),
            'all':( .5  , .5 ,5.   ),
            'ab': ( .5  , .5 ,1.e-9),
            'bg': (1.e-9, .5 ,5.   ),
            'ag': ( .5  ,1.e9,5.   ) }

for pkey, objval in paramlist.iteritems():
    for idrad in allradlist:
       nessfile = 'vesselness%s.%3.1f.nii.gz' % (pkey,idrad) 
       nesscmd = '/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artliver.nii.gz %s  1 1 %f %f %f %f' % (nessfile,objval[0],objval[1],objval[2],idrad)
       print(nesscmd)
       #if not os.path.isfile(nessfile):
       os.system(nesscmd)
       otsucmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter %s otsu%s.%3.1f.nii.gz 1  0' % (nessfile,pkey,idrad)
       print(otsucmd)
       os.system(otsucmd)

for pkey, objval in paramlist.iteritems():
    for (idlist,radlist) in enumerate(radlistlist):
      nessniilist =  ' '.join(['vesselness%s.%3.1f.nii.gz'%(pkey,idrad) for idrad in radlist])
      maxcmd     = 'c3d -verbose %s  -accum -max -endaccum  -o vesselmax%s%d.nii.gz' %(nessniilist ,pkey,idlist)
      print(maxcmd)
      os.system(maxcmd)
      otsumaxcmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter vesselmax%s%d.nii.gz otsumax%s%d.nii.gz 1  0' % (pkey,idlist,pkey,idlist)
      print(otsumaxcmd)
      os.system(otsumaxcmd)
      otsuniilist =  ' '.join(['otsu%s.%3.1f.nii.gz'%(pkey,idrad) for idrad in radlist])
      addbincmd = 'c3d -verbose %s -accum -add -endaccum -binarize -dilate 1 3x3x1vox -erode 1 3x3x1vox -o vessel%s%d.nii.gz' %(otsuniilist,pkey,idlist)
      print(addbincmd)
      os.system(addbincmd)
      overbincmd ='c3d -verbose manualregion.nii.gz  vessel%s%d.nii.gz -overlap 1 > overlap%s%02d.txt' %(pkey,idlist,pkey,idlist)
      overmaxcmd ='c3d -verbose manualregion.nii.gz otsumax%s%d.nii.gz -overlap 1 > overmax%s%02d.txt' %(pkey,idlist,pkey,idlist)
      print(overbincmd )
      print(overmaxcmd )
      os.system(overbincmd )
      os.system(overmaxcmd )


