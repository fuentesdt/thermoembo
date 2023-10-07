import subprocess
import os

allradlist = [1,2,3]
paramlist ={'dfl0':  (.5, .5 ,5.,1,1.,0 ),
            'dfl1':  (.5, .5 ,5.,1,1.,10),
            'dfl2':  (.5, .5 ,5.,1,1.,50),
            'dfl3':  (.5, .5 ,5.,1,.1,10),
            'dfl4':  (.5, .5 ,5.,1,.1,50),
            'dfl5':  (.5, .5 ,5.,1,10.,10),
            'dfl7':  (.5, .5 ,5.,1,10.,50),
#            'alp':  (.5,1.e9,0.,1,,),
#            'bet':  (0., .5 ,0.,1,,),
#            'gam':  (0.,1.e9,1.,1,,),
#            'all':  (.5, .5 ,1.,1,,),
#            'ab':   (.5, .5 ,0.,1,,),
#            'bg':   (0., .5 ,1.,1,,),
#            'ag':   (.5,1.e9,1.,1,,),
#            'nsalp':(.5   ,1.e9,0.,0),
#            'nsbet':(0., .5 ,0.,0),
#            'nsgam':(0.,1.e9,5.   ,0),
#            'nsall':( .5  , .5 ,5.   ,0),
#            'nsab': ( .5  , .5 ,0.   ,0),
#            'nsbg': (0.   , .5 ,5.   ,0),
#            'nsag': ( .5  ,1.e9,5.   ,0),
}
pixelsize = 1.0
pixelsize = 0.818359
maxhessdictionary = {}

# need hessian magnitude to scale gamma parameter for vesselness feature
# for idrad in allradlist:
#    hesscmd ='c3d -verbose artregion.nii.gz -hesseig %f -oo eig%d%%02d.nii.gz -foreach -dup -times -endfor -accum -add -endaccum -sqrt -o hessmag%d.nii.gz '%(idrad*pixelsize,idrad,idrad)
#    print(hesscmd)
#    #os.system(hesscmd)
#    getHeaderCmd = 'c3d hessmag%d.nii.gz liverregion.nii.gz  -lstat  ' % (idrad)
#    print (getHeaderCmd)
#    os.system( getHeaderCmd )
#    headerProcess = subprocess.Popen(getHeaderCmd ,shell=True,stdout=subprocess.PIPE )
#    while ( headerProcess.poll() == None ):
#       pass
#    rawlstatheader = filter(len,headerProcess.stdout.readline().strip('\n').split(" "))
#    rawlstatinfo = [filter(len,lines.strip('\n').split(" ")) for lines in headerProcess.stdout.readlines()]
#    labeldictionary =  dict([(int(line[0]),dict(zip(rawlstatheader[1:-1],map(float,line[1:-3])))) for line in rawlstatinfo ])
#    maxhessdictionary[idrad]= labeldictionary
# 
# print(maxhessdictionary)
   
for pkey, objval in paramlist.iteritems():
    denoisecmd = 'docker run --entrypoint=/opt/mitk/ImageDenoising/TotalVariationDenoisingImage  --rm -it --user $(id -u):$(id -g) -v /rsrch3/ip/dtfuentes/github/thermoembo/oncopig/zpaf22s016/:/data/  -v /rsrch3/ip/dtfuentes/github/thermoembo/oncopig/zpaf22s016:/out iptools:latest /data/artregion.nii.gz /out/artdenoise%3.1f%d.nii.gz %3.1f %d' % (objval[4],objval[5],objval[4],objval[5])
    print(denoisecmd)
    os.system(denoisecmd)
    for idrad in allradlist:
       nessfile = 'vesselness%s.%d.nii.gz' % (pkey,idrad) 
       #if pkey == 'dfl':
       nesscmd = '/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artdenoise%3.1f%d.nii.gz %s  1 1 %f %f %f %f %d' % (objval[4],objval[5],nessfile,objval[0],objval[1],objval[2],idrad*pixelsize,objval[3])
       #else:
       #   nesscmd = '/rsrch1/ip/dtfuentes/github/ExLib/Vesselness/HessianToObjectnessMeasureImageFilter artregion.nii.gz %s  1 1 %f %f %f %f %d' % (nessfile,objval[0],objval[1],objval[2]*maxhessdictionary[idrad][1]['Max']/2.,idrad*pixelsize,objval[3])
       print(nesscmd)
       #if not os.path.isfile(nessfile):
       os.system(nesscmd)
       otsucmd = 'c3d -verbose liverregion.nii.gz %s -times -o otsu%s.%d.nii.gz; /rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter otsu%s.%d.nii.gz  otsu%s.%d.nii.gz 1  0' % (nessfile,pkey,idrad,pkey,idrad,pkey,idrad)
       print(otsucmd)
       os.system(otsucmd)

radlistlist = [ allradlist[0:3], allradlist, allradlist[1:4]]
radlistlist = [ allradlist[0:3]]
for pkey, objval in paramlist.iteritems():
    for (idlist,radlist) in enumerate(radlistlist):
      nessniilist =  ' '.join(['vesselness%s.%d.nii.gz'%(pkey,idrad) for idrad in radlist])
      maxcmd     = 'c3d -verbose %s  -accum -max -endaccum liverregion.nii.gz -times -o vesselmax%s%d.nii.gz' %(nessniilist ,pkey,idlist)
      print(maxcmd)
      os.system(maxcmd)
      otsumaxcmd = '/rsrch1/ip/dtfuentes/github/ExLib/OtsuFilter/OtsuThresholdImageFilter vesselmax%s%d.nii.gz otsumax%s%d.nii.gz 1  0; c3d -verbose otsumax%s%d.nii.gz -dilate 1 3x3x1vox -erode 1 3x3x1vox -o otsumax%s%d.nii.gz' % (pkey,idlist,pkey,idlist,pkey,idlist,pkey,idlist)
      print(otsumaxcmd)
      os.system(otsumaxcmd)
      otsuniilist =  ' '.join(['otsu%s.%d.nii.gz'%(pkey,idrad) for idrad in radlist])
      addbincmd = 'c3d -verbose %s -accum -add -endaccum -binarize -dilate 1 3x3x1vox -erode 1 3x3x1vox -o vessel%s%d.nii.gz' %(otsuniilist,pkey,idlist)
      print(addbincmd)
      os.system(addbincmd)
      overbincmd ='c3d -verbose manualregion.nii.gz  vessel%s%d.nii.gz -overlap 1 > overlap%s%02d.txt' %(pkey,idlist,pkey,idlist)
      overmaxcmd ='c3d -verbose manualregion.nii.gz otsumax%s%d.nii.gz -overlap 1 > overmax%s%02d.txt' %(pkey,idlist,pkey,idlist)
      print(overbincmd )
      print(overmaxcmd )
      os.system(overbincmd )
      os.system(overmaxcmd )


