
    %% Load paths.
    if ~isdeployed
%      addpath('./common');
      addpath('./ooMFGRE');
%      addpath('./nifti');
    end

% /FUS4/data2/CJM/gitMATLAB/projects/mfgreScripts/
% /FUS4/data2/sfholtz/bitBucket/cressman/MGE_userguide_CJM.m
% /FUS4/data2/sfholtz/bitBucket/cressman/Mar2018Pub/SaveMGE_script.m



outputPathList = {'/FUS4/data2/sfholtz/Cressman/Mar2018Pub/EXP38_kidney1_03022017/MGE_obj_ARMA90_Model1.mat','/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney1Left_04202017_Exp42/MGE_obj_ARMA90_Model1.mat', '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/EXP38_kidney2_03022017/MGE_obj_ARMA90_Model1.mat','/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney2_01242017_Exp31/MGE_obj_ARMA90_Model1.mat', '/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney1_01242017_Exp31/MGE_obj_ARMA90_Model1.mat','/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney2Right_04202017_Exp42/MGE_obj_ARMA90_Model1.mat'}


outputPathList = {'/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney1Left_04202017_Exp42/newDeltaT.mat'}


for jjj = 1 : length(outputPathList )
    outputPath = outputPathList{jjj}
    fullpathsplit =  strsplit(outputPath  ,'/');
    outputdir = fullpathsplit{ length(fullpathsplit)-1}
    mkdir(outputdir)

    fulldata = load(outputPath);
    temperaturedata = fulldata.img.deltaT;
    
    
    header.PixelSpacing   = fulldata.img.hdrEx.PixelSpacing;
    header.SliceThickness = fulldata.img.sliceThickness ;
    
    nsteps = size(temperaturedata ,4);
    
    for iii = 1:nsteps
        fileout = sprintf('%s/temperature.%04d.vtk',outputdir ,iii)
        saveVTKFile(fileout ,'temperature', temperaturedata(:,:,:,iii), header)
    end

end





