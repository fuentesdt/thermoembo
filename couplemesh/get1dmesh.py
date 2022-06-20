
import vtk
# echo vtk version info
print("using vtk version", vtk.vtkVersion.GetVTKVersion())
import os

# load data
vtkReader = vtk.vtkPolyDataReader()
vtkReader.SetFileName( "CenterlineComputationModel.vtk" )

# merge duplicates data
cleaner = vtk.vtkCleanPolyData()
cleaner.SetInputConnection(vtkReader.GetOutputPort())
cleaner.Update()



DataVTK = cleaner.GetOutput()
points =DataVTK.GetPoints()
npoints =DataVTK.GetNumberOfPoints()
print(points.GetNumberOfPoints())

inVerts = DataVTK.GetVerts()
inLines = DataVTK.GetLines()
inPolys = DataVTK.GetPolys()
inStrips = DataVTK.GetStrips()

idList = vtk.vtkIdList()
cellNum = inLines.GetNumberOfCells()
cellLocation = 0
for iii in range(cellNum ):
  inLines.GetCell(cellLocation,idList )
  numnode = idList.GetNumberOfIds()
  # https://public.kitware.com/pipermail/vtkusers/2013-January/078002.html
  cellLocation += 1 + numnode ;
  print("numnode", numnode )
  #for jjj in range(numnode ):
  #  print(idList.GetId(jjj))



 
