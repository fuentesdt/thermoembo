# uses pyradiomics to generate feature maps and prints a histogram of the data values

# Requires pyradiomics version 2.0 or greater

# example usage: python createhistogram.py --file hackhistogram.csv 

from radiomics import featureextractor
import matplotlib.pyplot as plt
import six
import sys, os, argparse, csv
import numpy as np

import SimpleITK as sitk


# Get params from cmd line
parser = argparse.ArgumentParser(description = "Plot histograms of radiomics feature" +
    "for multiple images. Optionally save feature maps of the form imageName_maskName_labelNumber_featureName.nii.gz. Feature maps are cropped and need to be put back into the image space with c3d")
parser.add_argument("--file", "-f", help = "CSV with filepaths and columns Image, Mask, Label",
      required = True)
#parser.add_argument("--feature", default = None,
#      help = "feature name: i.e. original_firstorder_Mean. See pyradiomics -h for help")
parser.add_argument("--writeImage", "-w", dest="writeImage", action="store_true",
      help = "Include flag to write out a feature map (NIFTI) for each feature for each row. Cropped to the specific ROI")
parser.add_argument("--directory", "-d", default = "",
      help = "Directory to place output files in")
parser.add_argument("--param", "-p", default = "",
      help = "Pyradiomics parameter file to use for feature extraction. Old parameter files will need the VoxelSettings category added in order to be compatible.")

args = parser.parse_args()

# Check that csv file exists
if not os.path.isfile(args.file):
  sys.exit("Error: csv file " + args.file + " not found")

#Read file
inputCsv = csv.DictReader(open(args.file))
if 'Image' not in inputCsv.fieldnames:
  sys.exit("CSV file does not contain 'Image' column")
elif 'Mask' not in inputCsv.fieldnames:
  sys.exit("CSV file does not contain 'Mask' column")
elif 'Label' not in inputCsv.fieldnames:
  sys.exit("CSV file does not contain 'Label' column")


for idrow,row in enumerate(inputCsv):
  imageName = row['Image']
  maskName  = row['Mask']
  label     = int(row['Label'])

  # compute feature maps
  print("Extracting features for:\n" + imageName + "\n" + maskName + "\n" +  "label =" + str(label))
  ## result = extractor.execute(imageName, maskName, label, voxelBased=True)

  itkimage = sitk.ReadImage(imageName)
  itkmask  = sitk.ReadImage(maskName)

  print("Creating histogram " )
  voxelArray = sitk.GetArrayFromImage(itkimage)
  maskArray  = sitk.GetArrayFromImage(itkmask)
  #voxelArray = voxelArray[~np.isnan(voxelArray)] #only non-nan voxels
  voxelArray = voxelArray[np.nonzero( maskArray  == label )] #only get labels
  #voxelArray = voxelArray[np.argwhere( voxelArray > .05  )] 
  voxelMean = voxelArray.mean()
  print("Mean: " + str(voxelMean))

  labeldict = {1:'TE',2:'control',3:'emulsion'}
  #Print histogram
  if ( idrow % 3 == 0  ) :
    imageid =  imageName.split('/')[-1].split('.')[0] 
    outname = imageName.split('/')[-2] + imageName.split('/')[-1].split('.')[0] + maskName.split('/')[-1].split('.')[0] + "idrow%d" % idrow 
    print outname 
    plt.figure()
    plt.xlabel('intensity')
    if imageid == 'anatomy':
      plt.xlim(-100,500)
    if imageid == 'hessobj':
      plt.xlim(-10,100)
    plt.ylabel('frequency')
    plt.title('histogram')
  # HACK
  n, bins, patched = plt.hist(voxelArray,50,density=True,alpha=0.75, label = labeldict[label])
  if ( idrow % 3 == 2  ) :
    plt.grid(True)
    plt.legend()
    plt.savefig(outname  + '.png', bbox_inches='tight')

  ##   if args.writeImage:
  ##     print('Writing Image ' + outname)
  ##     sitk.WriteImage(val, outname)
 

# Check that feature name is valid


# featureExtractor.enableFeaturesByName()
## key is the feature class name, value is a list of enabled feature names
