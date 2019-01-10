
    %% Load paths.
    if ~isdeployed
%      addpath('./common');
      addpath('./ooMFGRE');
%      addpath('./nifti');
    end

% /FUS4/data2/CJM/gitMATLAB/projects/mfgreScripts/
% /FUS4/data2/sfholtz/bitBucket/cressman/MGE_userguide_CJM.m
% /FUS4/data2/sfholtz/bitBucket/cressman/Mar2018Pub/SaveMGE_script.m

outputPath = '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney1_01122017_Exp29/MGE_obj.mat'; % Path to write MGE object to disk
outputPath = '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney1_01242017_Exp31/MGE_obj.mat'; % Path to write MGE object to disk
outputPath = '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney2_01242017_Exp31/MGE_obj.mat'; % Path to write MGE object to disk
outputPath = '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/EXP38_kidney1_03022017/MGE_obj.mat'; % Path to write MGE object to disk
outputPath = '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney1Left_04202017_Exp42/MGE_obj.mat'; % Path to write MGE object to disk
outputPath = '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney2Right_04202017_Exp42/MGE_obj.mat'; % Path to write MGE object to disk
outputPath = '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/EXP38_kidney2_03022017/MGE_obj.mat'; % Path to write MGE object to disk

outputPath = '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney2_01242017_Exp31/MGE_obj_ARMA90_Model1.mat';


xxx = load(outputPath)


yyy = xxx.img.deltaT;


header.PixelSpacing   = xxx.img.hdrEx.PixelSpacing;
header.SliceThickness = xxx.img.sliceThickness ;

nsteps = size(yyy,4);

for iii = 1:nsteps
    fileout = sprintf('temperature.%04d.vtk',iii)
    saveVTKFile(fileout ,'temperature', yyy(:,:,:,iii), header)
end




