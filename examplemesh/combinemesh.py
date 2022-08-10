# load centerline file and merge with background sphere mesh
import os
import numpy as np
import nibabel as nib  
import vtk
import ctypes 

# c3d temperature.0050.vtk -replace NaN 0 maxtemp.vtk

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--file_name",
                  action="store", dest="file_name", default="meshSphere.e",
                  help="getting nodes from this file", metavar = "FILE")
parser.add_option( "--output",
                  action="store", dest="output", default="mytetmesh",
                  help="converting/this/file to mesh", metavar = "FILE")
(options, args) = parser.parse_args()
if (options.file_name):
  
  ## load nodes
  vtkExodusIIReader = vtk.vtkExodusIIReader()
  vtkExodusIIReader.SetFileName( options.file_name )
  vtkExodusIIReader.Update()	  
  exomesh =  vtkExodusIIReader.GetOutput() 
  exomesh.IsA("vtkMultiBlockDataSet")
  myiter = exomesh.NewIterator()
  myiter.UnRegister(None)
  myiter.InitTraversal()
  curInput = myiter.GetCurrentDataObject()
  mynodes = curInput.GetPoints()
    
  # loop over blocks...
  #while not myiter.IsDoneWithTraversal():
  #    curInput = myiter.GetCurrentDataObject()

  #exomesh.GetNumberOfElements()

  # FIXME notice that order of operations is IMPORTANT
  # FIXME translation followed by rotation will give different results
  # FIXME than rotation followed by translation
  # FIXME Translate -> RotateZ -> RotateY -> RotateX -> Scale seems to be the order of paraview
  AffineTransform = vtk.vtkTransform()
  # should be in meters
  convertmm = 0.001
  AffineTransform.Translate(0,0,0 )
  AffineTransform.RotateZ( 0 )
  AffineTransform.RotateY( 0 )
  AffineTransform.RotateX( 0 )
  AffineTransform.Scale([convertmm ,convertmm ,convertmm ])

  ## load vessel
  vesselreader = vtk.vtkDataSetReader();
  vesselreader.SetFileName("Centerlinemodel.vtk");
  vesselreader.Update();
  vtkTransformFEMMesh = vtk.vtkTransformFilter()
  vtkTransformFEMMesh.SetTransform( AffineTransform )
  vtkTransformFEMMesh.SetInputData( vesselreader.GetOutput() )
  vtkTransformFEMMesh.Update()
  ## remove nodes that are close/duplicates
  vtkClean = vtk.vtkCleanPolyData()
  vtkClean.SetTolerance(0.05)
  vtkClean.SetInputData( vtkTransformFEMMesh.GetOutput() )
  vtkClean.Update( )
  VesselDataVMTK = vtkClean.GetOutput();
  ## convert from polydata to line elements
  points  = VesselDataVMTK.GetPoints()
  inLines = VesselDataVMTK.GetLines()
  inField = VesselDataVMTK.GetPointData()
  cellLocation = 0
  # use set to remove duplicate
  cellSet = set()
  idList = vtk.vtkIdList()
  cellNum = inLines.GetNumberOfCells()
  for iii in range(cellNum ):
    inLines.GetCell(cellLocation,idList )
    numnode = idList.GetNumberOfIds()
    # https://public.kitware.com/pipermail/vtkusers/2013-January/078002.html
    cellLocation += 1 + numnode ;
    print("numnode", numnode )
    for jjj in range(numnode-1):
      cellSet.add((idList.GetId(jjj),idList.GetId(jjj+1) ))
      #print(idList.GetId(jjj))
  
  # use set to remove duplicate
  lineElements= vtk.vtkCellArray()
  for idline in cellSet:
      myLineElement = vtk.vtkIdList()
      #print(idline)
      myLineElement.InsertNextId(idline[0])
      myLineElement.InsertNextId(idline[1])
      lineElements.InsertNextCell(myLineElement )
  
  # set polydata
  VesselData = vtk.vtkPolyData()
  VesselData.SetPoints(points)
  VesselData.GetPointData().AddArray(inField.GetArray(0))
  VesselData.SetVerts( lineElements )

  vesselwriter = vtk.vtkDataSetWriter();
  vesselwriter.SetFileName("Centerlinemodeltransform.vtk" )
  vesselwriter.SetInputData(VesselData )
  vesselwriter.Update()

  nodes=[]
  for iii in range(mynodes.GetNumberOfPoints() ):
    nodes.append(mynodes.GetPoint(iii))
  
  coarsenode=[]
  for iii in range(mynodes.GetNumberOfPoints() ):
    coarsenode.append(mynodes.GetPoint(iii))

  # add surfacepoints
  for ipoint in range( VesselData.GetPoints().GetNumberOfPoints() ):
      CurrentPoint = VesselData.GetPoint(ipoint)
      newnode = ( CurrentPoint[0] ,
                  CurrentPoint[1] ,
                  CurrentPoint[2] )
      coarsenode.append(newnode)

  ## 
  ## distance = numpyimage[:,:,:,:,3]
  ## distancemax = np.max(distance[:])
  ## distancemin = max(np.min(distance[:]),0)
  ## 
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
    #f.write("%12.5e\n" % ( max(pixelsize[0],(1. - np.exp(-5*node[3]/distancemax)) * 32.0 * pixelsize[0] ) ) )
    f.write("%12.5e\n" % 1. )
  f.close()
  
  f = open("%s.2.b.mtr" % options.output, "w")
  f.write("%d 1 \n" % len(nodes))
  for iii, node in enumerate(nodes):
    #f.write("%12.5e\n" % ( max(0.5*pixelsize[0],(1. - np.exp(-4*node[3]/distancemax)) * 32.0 * pixelsize[0] ) ) )
    f.write("%12.5e\n" % 1. )
  f.close()
  
  f = open("%s.3.b.mtr" % options.output, "w")
  f.write("%d 1 \n" % len(nodes))
  for iii, node in enumerate(nodes):
    #f.write("%12.5e\n" % ( max(0.5*pixelsize[0],(1. - np.exp(-3*node[3]/distancemax)) * 32.0 * pixelsize[0] ) ) )
    f.write("%12.5e\n" % 1. )
  f.close()
  
  
  # run adaptive mesh generation
  backgroundmeshcmd = "tetgen -k %sbg.node;ln -sf %sbg.1.node  %s.1.b.node ; ln -sf %sbg.1.ele   %s.1.b.ele ;ln -sf %sbg.1.node  %s.2.b.node ; ln -sf %sbg.1.ele   %s.2.b.ele ;ln -sf %sbg.1.node  %s.3.b.node ; ln -sf %sbg.1.ele   %s.3.b.ele " % (options.output, options.output,options.output, options.output,options.output, options.output,options.output, options.output,options.output, options.output,options.output, options.output,options.output)
  initializecoarsemeshcmd = "tetgen -k %s.node"       % options.output
  #adaptivemeshcmd = "tetgen -rmqk %s.1.node;tetgen -rmqk %s.2.node; tetgen -rmqk %s.3.node" % (options.output, options.output, options.output)
  adaptivemeshcmd = "tetgen -rmqk %s.1.node" % options.output
                    
  print(backgroundmeshcmd )
  print(initializecoarsemeshcmd )
  print(adaptivemeshcmd )
  os.system(backgroundmeshcmd )
  os.system(initializecoarsemeshcmd )
  #os.system(adaptivemeshcmd )

  
  # load vtk data
  #for iii in [1,2,3,4]:
  for iii in [1]:
     vtkReader = vtk.vtkUnstructuredGridReader()
     vtkReader.SetFileName( "%s.%d.vtk" % (options.output,iii) )
     vtkReader.Update()
     vtkMesh  = vtkReader.GetOutput()
     print('setup metadata' ,iii)

     mytol = 1.e-5
     vesselmap = {}
     # FIXME - HACK - need to quick sort then find 
     for ipoint in range( vtkMesh  .GetPoints().GetNumberOfPoints() ):
         CurrentPoint = vtkMesh.GetPoint(ipoint)
         newnode = ( CurrentPoint[0] ,
                     CurrentPoint[1] ,
                     CurrentPoint[2] )
         for jpoint in range( VesselData.GetPoints().GetNumberOfPoints() ):
             vPoint = VesselData.GetPoint(jpoint)
             vnode = ( vPoint[0] ,
                         vPoint[1] ,
                         vPoint[2] )
             if( abs(vPoint[0] - CurrentPoint[0]) < mytol   
             and abs(vPoint[1] - CurrentPoint[1]) < mytol   
             and abs(vPoint[2] - CurrentPoint[2]) < mytol ):
               print(ipoint,jpoint,newnode,vnode)
               vesselmap[jpoint] = ipoint
     for idset, idline in enumerate (cellSet):
         print(idset,idline ,vesselmap[idline[0]],vesselmap[idline[1]])

     ## vtkNew<vtkDummyController> controller;
     ## controller->Initialize(&argc, &argv, 1);
     ## vtkMultiProcessController::SetGlobalController(controller.Get());
     ## HACK for parallel write
     ## https://www.paraview.org/Bug/view.php?id=15813
     controller =vtk.vtkDummyController()
     vtk.vtkMultiProcessController.SetGlobalController(controller)

     print('convert to exodus' ,iii)
     # convert to exodus 
     myexowriter = vtk.vtkExodusIIWriter()
     myexowriter.DebugOn()
     myexowriter.SetFileName("%s.%d.exo" % (options.output,iii))
     myexowriter.SetInputData( vtkMesh )
     print('writing 1' ,"%s.%d.exo" % (options.output,iii))
     myexowriter.Update()
     ## FIXME - need to setup node sets
     ## ## # create a 3-int array
     ## ## pynodelist = [1,2,3,4,5]
     ## ## cnodelist = (ctypes.c_int*5)(*pynodelist )
     ## ## 
     ## ## cnodesetid =ctypes.c_int(4)
     ## ## cnodesetsize =ctypes.c_int(5)

     ## testem = myexowriter.GetModelMetadata();
     ## testem.SetTitle("dftest")
     ## testem.SetNumberOfNodeSets(1)
     ## print testem.GetDimension();
     ## print testem.GetNumberOfNodeSets();
     ## testem.SetNodeSetSize(ctypes.byref(cnodesetsize))
     ## testem.SetNodeSetIds(ctypes.pointer(cnodesetid))
     ## testem.SetNodeSetNodeIdList(cnodelist )
     ## print  myexowriter
     ## print 'writing 1' ,"%s.%d.exo" % (options.output,iii)
     ## myexowriter.Update()




else:
  parser.print_help()
  print(options)
## with open("isosurface.1.node") as myfile:
##     lines = [list(filter(len,line.strip().split(' '))) for line in myfile]
## numnodes = int(lines[0][0] )

