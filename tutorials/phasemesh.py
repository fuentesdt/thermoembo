import os
import numpy as np
import nibabel as nib  
import vtk


# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--file_name",
                  action="store", dest="file_name", default="../temperaturedata/Kidney1Left_04202017_Exp42/vessel.nii.gz",
                  help="converting/this/file to mesh", metavar = "FILE")
parser.add_option( "--output",
                  action="store", dest="output", default="mytetmesh",
                  help="converting/this/file to mesh", metavar = "FILE")
(options, args) = parser.parse_args()
if (options.file_name):
  c3dcmd = "c3d -verbose %s  -o %simage.vtk  -sdt -as A -cmp -push A -omc %ssetup.nii.gz" % (options.file_name,options.output,options.output)
  print c3dcmd 
  os.system( c3dcmd )
  
  imagedata = nib.load("%ssetup.nii.gz" % options.output)
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
  
  f = open("%sbg.node" % options.output, "w")
  f.write("%d 3 0 0\n" % len(nodes))
  for iii, node in enumerate(nodes):
    f.write("%d %f %f %f\n" % (iii, node[0],node[1],node[2]))
  f.close()
  
  f = open("%s.node" % options.output, "w")
  f.write("%d 3 0 0\n" % len(coarsenode))
  for iii, node in enumerate(coarsenode):
    f.write("%d %f %f %f\n" % (iii, node[0],node[1],node[2]))
  f.close()
  
  
  f = open("%s.1.b.mtr" % options.output, "w")
  f.write("%d 1 \n" % len(nodes))
  for iii, node in enumerate(nodes):
    f.write("%f\n" % ( max(zooms[0],node[3]/distancemax * 25.0 * zooms[0]) ) )
    #f.write("1.0\n"                                                    )
  f.close()
  
  # run adaptive mesh generation
  backgroundmeshcmd = "tetgen -k %sbg.node;cp %sbg.1.node  %s.1.b.node ; cp %sbg.1.ele   %s.1.b.ele " % (options.output, options.output,options.output, options.output,options.output)
  initializecoarsemeshcmd = "tetgen -k %s.node"       % options.output
  adaptivemeshcmd = "tetgen -Vrmqk %s.1.node" % options.output
  print backgroundmeshcmd 
  print initializecoarsemeshcmd 
  print adaptivemeshcmd 
  #os.system(backgroundmeshcmd )
  #os.system(initializecoarsemeshcmd )
  #os.system(adaptivemeshcmd )

  # load vtk data
  vtkReader = vtk.vtkUnstructuredGridReader()
  vtkReader.SetFileName( "%s.2.vtk" % options.output )
  vtkReader.Update()
  xxx = vtkReader.GetOutput() 
  yyy = vtkReader.GetOutputPort() 

  # convert to exodus 
  vtkExodusIIWriter = vtk.vtkExodusIIWriter()
  vtkExodusIIWriter.SetFileName("%s.2.e" % options.output)
  vtkExodusIIWriter.SetInputConnection( vtkReader.GetOutputPort()  )
  vtkExodusIIWriter.WriteAllTimeStepsOn()
  vtkExodusIIWriter.Write()

else:
  parser.print_help()
  print options
## with open("isosurface.1.node") as myfile:
##     lines = [list(filter(len,line.strip().split(' '))) for line in myfile]
## numnodes = int(lines[0][0] )

