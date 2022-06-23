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
cellNum = inLines.GetNumberOfCells()
cellLocation = 0
# use set to remove duplicate
cellSet = set()
idList = vtk.vtkIdList()
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
    print(idline)
    myLineElement.InsertNextId(idline[0])
    myLineElement.InsertNextId(idline[1])
    lineElements.InsertNextCell(myLineElement )

# set polydata
polydata = vtk.vtkPolyData()
polydata.SetPoints(points)
polydata.SetVerts( lineElements )

OutputFileName="testline.vtk"
# write to file
polydatawriter = vtk.vtkDataSetWriter()
polydatawriter.SetFileName(OutputFileName)
polydatawriter.SetInputData(polydata)
polydatawriter.Update()
