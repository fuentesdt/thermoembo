import os
import vtk
import vtk.util.numpy_support as vtkNumPy 
import numpy as np
import nibabel as nib  


with open("isosurface.1.node") as myfile:
    lines = [list(filter(len,line.strip().split(' '))) for line in myfile]

numnodes = int(lines[0][0] )

## vtkReader = vtk.vtkDataSetReader() 
## vtkReader.SetFileName( "signed.vtk" ) 
## vtkReader.Update()
## imageDataVTK = vtkReader.GetOutput()
## dimensions = imageDataVTK.GetDimensions()
## spacing = imageDataVTK.GetSpacing()
## origin  = imageDataVTK.GetOrigin()
## print spacing, origin, dimensions
## #fem.SetImagingDimensions( dimensions ,origin,spacing) 
## 
## image_point_data = imageDataVTK.GetPointData() 
## image_data       = vtkNumPy.vtk_to_numpy( image_point_data.GetArray(0) ) 
import nibabel as nib  
imagedata = nib.load("test.nii.gz" )
numpyimage= imagedata.get_data()

nodes=[]
## for iii in  range(numnodes):
##    nodes.append( (float(lines[iii+1][1]),float(lines[iii+1][2]),float(lines[iii+1][3]), .1) )

#image data
dim = numpyimage.shape

for iii in range(dim[0]):
  for jjj in range(dim[1]):
     for kkk in range(dim[2]):
         newnode = (numpyimage[iii,jjj,kkk,0,0],
                    numpyimage[iii,jjj,kkk,0,1],
                    numpyimage[iii,jjj,kkk,0,2],
                    numpyimage[iii,jjj,kkk,0,3])
         nodes.append(newnode)

f = open("mytest.node", "w")
f.write("%d 3 0 0\n" % len(nodes))
for iii, node in enumerate(nodes):
  f.write("%d %f %f %f\n" % (iii, node[0],node[1],node[2]))
f.close()

f = open("mytest.mrt", "w")
f.write("%d 1 \n" % len(nodes))
for iii, node in enumerate(nodes):
  f.write("%f\n" % (node[3]))
f.close()

os.system("tetgen -k mytest.node")

