
# /opt/apps/slicer/Slicer-5.0.2-linux-amd64/Slicer --launcher-verbose --python-script "/rsrch3/ip/dtfuentes/github/thermoembo/oncopig/zpaf23s023/autocenterline.py" --no-splash --no-main-window
loadedVolumeNode = slicer.util.loadVolume("/rsrch3/ip/dtfuentes/github/thermoembo/oncopig/zpaf23s023/art.nii.gz")


# Create segmentation from a NIFTI + color table file
#colorNode = slicer.util.loadColorTable('c:/tmp/tmp/Segmentation-label_ColorTable.ctbl')
slicer.util.loadSegmentation("/rsrch3/ip/dtfuentes/github/thermoembo/oncopig/zpaf23s023/ven.manual.nii.gz")


segmentationNode = slicer.mrmlScene.GetNodeByID("vtkMRMLSegmentationNode1")
referenceVolumeNode = slicer.mrmlScene.GetNodeByID("vtkMRMLScalarVolumeNode1")
labelmapVolumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode")
# #segmentationNode = slicer.util.getNode("Segmentation")
# segmentNames = ["Segment_1", "Segment_2"]
# segmentIds = vtk.vtkStringArray()
# for segmentName in segmentNames:
#   segmentId = segmentationNode.GetSegmentation().GetSegmentIdBySegmentName(segmentName)
#   segmentIds.InsertNextValue(segmentId)
# slicer.vtkSlicerSegmentationsModuleLogic.ExportSegmentsToLabelmapNode(segmentationNode, segmentIds, labelmapVolumeNode, referenceVolumeNode)

#labelmapVolumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode")
#slicer.modules.segmentations.logic().ExportAllSegmentsToLabelmapNode(segmentationNode, labelmapVolumeNode, slicer.vtkSegmentation.EXTENT_REFERENCE_GEOMETRY)


shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
exportFolderItemId = shNode.CreateFolderItem(shNode.GetSceneItemID(), "Segments")
slicer.modules.segmentations.logic().ExportAllSegmentsToModels(segmentationNode, exportFolderItemId)

