#### track manual updated nifti label files with git lfs
https://git-lfs.github.com/

#### use slicer to extract centerline

ignore file orientation and/or export the slicer import images to register the centerline to the image data

./ignoreorientation.PNG

the vessel tips are typically too small and cause connection errors... dilarte the vessel data to have the algorithm find the center of the dilated vessel 

./keystepdilation.PNG
