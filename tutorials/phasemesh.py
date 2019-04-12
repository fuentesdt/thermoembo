import os
import numpy as np
import nibabel as nib  


## with open("isosurface.1.node") as myfile:
##     lines = [list(filter(len,line.strip().split(' '))) for line in myfile]
## numnodes = int(lines[0][0] )

## c3d -verbose ../temperaturedata/Kidney1Left_04202017_Exp42/vesselregion.vtk  -sdt -as A -cmp -push A -omc test.nii.gz
import nibabel as nib  
imagedata = nib.load("test.nii.gz" )
zooms = imagedata.header.get_zooms()
numpyimage= imagedata.get_data()

## for iii in  range(numnodes):
##    nodes.append( (float(lines[iii+1][1]),float(lines[iii+1][2]),float(lines[iii+1][3]), .1) )

#image data
dim = numpyimage.shape

nodes=[]
for iii in range(dim[0]):
  for jjj in range(dim[1]):
     for kkk in range(dim[2]):
         newnode = (numpyimage[iii,jjj,kkk,0,0],
                    numpyimage[iii,jjj,kkk,0,1],
                    numpyimage[iii,jjj,kkk,0,2],
                    numpyimage[iii,jjj,kkk,0,3])
         nodes.append(newnode)

coarsenode=[]
for iii in range(0,dim[0],10):
  for jjj in range(0,dim[1],10):
     for kkk in range(0,dim[2],2):
         newnode = (numpyimage[iii,jjj,kkk,0,0],
                    numpyimage[iii,jjj,kkk,0,1],
                    numpyimage[iii,jjj,kkk,0,2],
                    numpyimage[iii,jjj,kkk,0,3])
         coarsenode.append(newnode)

distance = numpyimage[:,:,:,:,3]
distancemax = np.max(distance[:])
distancemin = max(np.min(distance[:]),0)

f = open("mytest.node", "w")
f.write("%d 3 0 0\n" % len(nodes))
for iii, node in enumerate(nodes):
  f.write("%d %f %f %f\n" % (iii, node[0],node[1],node[2]))
f.close()

f = open("coarse.node", "w")
f.write("%d 3 0 0\n" % len(coarsenode))
for iii, node in enumerate(coarsenode):
  f.write("%d %f %f %f\n" % (iii, node[0],node[1],node[2]))
f.close()


f = open("coarse.1.b.mtr", "w")
f.write("%d 1 \n" % len(nodes))
for iii, node in enumerate(nodes):
  f.write("%f\n" % ( max(zooms[0],node[3]/distancemax * 20.0 * zooms[0]) ) )
  #f.write("1.0\n"                                                    )
f.close()

os.system("tetgen -k mytest.node")
os.system("cp mytest.1.node  coarse.1.b.node")
os.system("cp mytest.1.ele   coarse.1.b.ele ")
os.system("tetgen -k coarse.node")
os.system("tetgen -Vrmqk coarse.1.node")

