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
  c3dcmd = "c3d -verbose %s -o %simage.vtk -as B -dilate 1 1x1x0vox -push B  -scale -1 -add -as C -sdt -as A -cmv -push A -push C -omc %ssetup.nii.gz" % (options.file_name,options.output,options.output)
  print c3dcmd 
  os.system( c3dcmd )
  
  imagedata = nib.load("%ssetup.nii.gz" % options.output)
  zooms = imagedata.header.get_zooms()
  numpyimage= imagedata.get_data()
  
  binaryReader = vtk.vtkDataSetReader()
  binaryReader.SetFileName( "%simage.vtk" % options.output )
  binaryReader.Update()

  #surfExtractor = vtk.vtkMarchingCubes()
  surfExtractor = vtk.vtkContourFilter()
  surfExtractor.SetInputData(binaryReader.GetOutput())
  surfExtractor.SetValue(0,.7)
  surfExtractor.Update()
  surface = surfExtractor.GetOutput() 

  stlwriter = vtk.vtkSTLWriter()
  stlwriter.SetInputData( surface )
  stlwriter.SetFileName( "%simage.stl" % options.output )
  stlwriter.SetFileTypeToASCII()
  #stlwriter.SetFileTypeToBinary()
  stlwriter.Write()


  ## for iii in  range(numnodes):
  ##    nodes.append( (float(lines[iii+1][1]),float(lines[iii+1][2]),float(lines[iii+1][3]), .1) )
  
  #image data
  dim = numpyimage.shape

  mmconversion = 1.e-3
  pixelsize = (mmconversion* zooms[0] ,mmconversion* zooms[1] ,mmconversion* zooms[2] )
  
  nodes=[]
  for iii in range(dim[0]):
    for jjj in range(dim[1]):
       for kkk in range(dim[2]):
           newnode = (pixelsize[0] * numpyimage[iii,jjj,kkk,0,0],
                      pixelsize[1] * numpyimage[iii,jjj,kkk,0,1],
                      pixelsize[2] * numpyimage[iii,jjj,kkk,0,2],
                                     numpyimage[iii,jjj,kkk,0,3],
                                     numpyimage[iii,jjj,kkk,0,4])
           nodes.append(newnode)
  
  coarsenode=[]
  for iii in range(0,dim[0],32):
    for jjj in range(0,dim[1],32):
       for kkk in range(0,dim[2],2):
           newnode = (pixelsize[0] * numpyimage[iii,jjj,kkk,0,0],
                      pixelsize[1] * numpyimage[iii,jjj,kkk,0,1],
                      pixelsize[2] * numpyimage[iii,jjj,kkk,0,2])
           coarsenode.append(newnode)

  # add surfacepoints
  for ipoint in range( surface.GetPoints().GetNumberOfPoints() ):
      CurrentPoint = surface.GetPoint(ipoint)
      newnode = (mmconversion * CurrentPoint[0] ,
                 mmconversion * CurrentPoint[1] ,
                 mmconversion * CurrentPoint[2] )
      coarsenode.append(newnode)

  
  distance = numpyimage[:,:,:,:,3]
  boundaryid = numpyimage[:,:,:,:,4]
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
    #f.write("%12.5e\n" % ( 0.5*pixelsize[0]*node[4] ) )
    f.write("%12.5e\n" % ( max(pixelsize[0],(1. - np.exp(-5*node[3]/distancemax)) * 32.0 * pixelsize[0] ) ) )
  f.close()
  
  f = open("%s.2.b.mtr" % options.output, "w")
  f.write("%d 1 \n" % len(nodes))
  for iii, node in enumerate(nodes):
    f.write("%12.5e\n" % ( max(0.5*pixelsize[0],(1. - np.exp(-4*node[3]/distancemax)) * 32.0 * pixelsize[0] ) ) )
  f.close()
  
  f = open("%s.3.b.mtr" % options.output, "w")
  f.write("%d 1 \n" % len(nodes))
  for iii, node in enumerate(nodes):
    f.write("%12.5e\n" % ( max(0.5*pixelsize[0],(1. - np.exp(-3*node[3]/distancemax)) * 32.0 * pixelsize[0] ) ) )
  f.close()
  
  
  # run adaptive mesh generation
  backgroundmeshcmd = "tetgen -k %sbg.node;ln -sf %sbg.1.node  %s.1.b.node ; ln -sf %sbg.1.ele   %s.1.b.ele ;ln -sf %sbg.1.node  %s.2.b.node ; ln -sf %sbg.1.ele   %s.2.b.ele ;ln -sf %sbg.1.node  %s.3.b.node ; ln -sf %sbg.1.ele   %s.3.b.ele " % (options.output, options.output,options.output, options.output,options.output, options.output,options.output, options.output,options.output, options.output,options.output, options.output,options.output)
  initializecoarsemeshcmd = "tetgen -k %s.node"       % options.output
  adaptivemeshcmd = "tetgen -rmqk %s.1.node;tetgen -rmqk %s.2.node; tetgen -rmqk %s.3.node" % (options.output, options.output, options.output)
                    
  print backgroundmeshcmd 
  print initializecoarsemeshcmd 
  print adaptivemeshcmd 
  os.system(backgroundmeshcmd )
  os.system(initializecoarsemeshcmd )
  os.system(adaptivemeshcmd )

  
  # load vtk data
  for iii in [2,3,4]:
     vtkReader = vtk.vtkUnstructuredGridReader()
     vtkReader.SetFileName( "%s.%d.vtk" % (options.output,iii) )
     vtkReader.Update()


     ## vtkNew<vtkDummyController> controller;
     ## controller->Initialize(&argc, &argv, 1);
     ## vtkMultiProcessController::SetGlobalController(controller.Get());
     ## HACK for parallel write
     ## https://www.paraview.org/Bug/view.php?id=15813
     controller =vtk.vtkDummyController()
     vtk.vtkMultiProcessController.SetGlobalController(controller)

     # convert to exodus 
     vtkExodusIIWriter = vtk.vtkExodusIIWriter()
     #vtkExodusIIWriter.DebugOn()
     vtkExodusIIWriter.SetFileName("%s.%d.exo" % (options.output,iii))
     vtkExodusIIWriter.SetInputData( vtkReader.GetOutput() )
     print vtkExodusIIWriter
     vtkExodusIIWriter.Update()

else:
  parser.print_help()
  print options
## with open("isosurface.1.node") as myfile:
##     lines = [list(filter(len,line.strip().split(' '))) for line in myfile]
## numnodes = int(lines[0][0] )

