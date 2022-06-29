#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkStructuredPoints.h>
#include <vtkDataSetReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPolyData.h>
#include <vtkExodusIIWriter.h>
#include <vtkModelMetadata.h>

int main(int argc, char **argv)
{
  vtkSmartPointer<vtkDataSetReader> vesselreader = vtkSmartPointer<vtkDataSetReader>::New();

  vtkSmartPointer<vtkUnstructuredGridReader>  vtkReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
  vtkReader->SetFileName( "mytetmesh.2.vtk"  );
  vtkReader->Update();
  //vtkSmartPointer<vtkUnstructuredGrid>  mymesh = vtkReader->GetOutput();

  //## vtkNew<vtkDummyController> controller;
  //## controller->Initialize(&argc, &argv, 1);
  //## vtkMultiProcessController::SetGlobalController(controller.Get());
  //## HACK for parallel write
  //## https://www.paraview.org/Bug/view.php?id=15813
  //controller =vtk.vtkDummyController()
  //vtk.vtkMultiProcessController.SetGlobalController(controller)

  //   print 'convert to exodus' ,iii
  //   # convert to exodus 
  vtkSmartPointer<vtkExodusIIWriter> myexowriter = vtkSmartPointer<vtkExodusIIWriter>::New();
     myexowriter->DebugOn();
     myexowriter->SetFileName("mytetmesh.2.exo" );
     myexowriter->SetInputConnection( vtkReader->GetOutputPort());
     myexowriter->Update();
     //print 'writing 1' ,"%s.%d.exo" % (options.output,iii)
  vtkSmartPointer<vtkModelMetadata>   testem = myexowriter->GetModelMetadata();
     testem->SetTitle("dftest");
     testem->SetNumberOfNodeSets(1);
   std::cout << testem->GetDimension() << std::endl;
   std::cout << testem->GetNumberOfNodeSets() << std::endl;
   int cnodesetsize = 5;
   int cnodesetid = 4;
   testem->SetNodeSetSize(&cnodesetsize );
   testem->SetNodeSetIds( &cnodesetid  );
   //  testem.SetNodeSetNodeIdList(cnodelist )
   //  print  myexowriter
   myexowriter->PrintSelf(std::cout,vtkIndent());
   //  print 'writing 1' ,"%s.%d.exo" % (options.output,iii)
     myexowriter->Update();
  return 0;
}
