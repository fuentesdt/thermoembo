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

# get data structures
DataVTK = cleaner.GetOutput()
points =DataVTK.GetPoints()
npoints =DataVTK.GetNumberOfPoints()
print(points.GetNumberOfPoints())

inVerts = DataVTK.GetVerts()
inLines = DataVTK.GetLines()
inPolys = DataVTK.GetPolys()
inStrips = DataVTK.GetStrips()

# convert from polydata to line elements
idList = vtk.vtkIdList()
lineElements= vtk.vtkCellArray()
cellNum = inLines.GetNumberOfCells()
cellLocation = 0
for iii in range(cellNum ):
  inLines.GetCell(cellLocation,idList )
  numnode = idList.GetNumberOfIds()
  # https://public.kitware.com/pipermail/vtkusers/2013-January/078002.html
  cellLocation += 1 + numnode ;
  print("numnode", numnode )
  for jjj in range(numnode-1):
    myLineElement = vtk.vtkIdList()
    myLineElement.InsertNextId(idList.GetId(jjj))
    myLineElement.InsertNextId(idList.GetId(jjj+1))
    lineElements.InsertNextCell(myLineElement )
    #print(idList.GetId(jjj))

# set polydata
polydata = vtk.vtkPolyData()
polydata.SetPoints(points)
polydata.SetVerts( lineElements )

# merge duplicate lines 
linecleaner = vtk.vtkCleanPolyData()
linecleaner.SetInputData(polydata)
linecleaner.Update()

OutputFileName="testline.vtk"
# write to file
polydatawriter = vtk.vtkDataSetWriter()
polydatawriter.SetFileName(OutputFileName)
polydatawriter.SetInputData(linecleaner.GetOutput())
polydatawriter.Update()

 
