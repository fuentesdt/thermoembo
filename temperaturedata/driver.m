% works with matlab 2018 

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

outputPathList = {'/mnt/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney1Left_04202017_Exp42/newDeltaT.mat', '/mnt/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney2Right_04202017_Exp42/newDeltaT.mat' , '/mnt/FUS4/data2/sfholtz/Cressman/Mar2018Pub/EXP38_kidney2_03022017/newDeltaT.mat' }

outputPathList = {'/mnt/FUS4/data2/sfholtz/Cressman/Mar2018Pub/Kidney1Left_04202017_Exp42/newDeltaT.mat'}


% https://www.tandfonline.com/doi/full/10.1080/02656736.2019.1635274
% MRTI measurements were performed on a 3.0 T MR scanner (750 W, GE Healthcare Technologies, Waukesha, WI, USA) using the body coil for signal excitation with a 32-channel phased-array cardiac coil. Temperature estimation was based on the water proton resonance frequency (PRF) shift of water [22], provided by a 2 D multiple fast gradient-recalled echoes (MFGRE) sequence [23,24] set to output real and imaginary data at each echo time. Parameters for 2 D MFGRE were as follows: time of repetition (TR)/First echo time (TE)/Last TE = 100/1.50/10 ms; 15 total echoes in 3 shots of 5 echoes; echo spacing = 0.604 ms; parallel imaging acceleration factor = 1.5; flip angle = 25°; acquisition matrix = 96 × 96; reconstructed matrix = 256 × 256; field of view = 180 × 180 mm2; number of slices = 5; slice spacing and thickness = 5 mm; receiver bandwidth =  ±83.3 kHz; spectral bandwidth = 165 kHz or 13.0 parts per million (ppm); the time per dynamic frame ranged from 20.46 s to 25.94 s in the five experiments.

echotime = 1.5e-3; % s
field    = 3.; % T
alpha    = .01; % ppm/degC
gamma    = 42.58; % MHz/T
echospacing=.604; %ms
numberecho = 15; 
% c3d Kidney1Left_04202017_Exp42/magnitude.0001.vtk -replace NaN 0 -o  Kidney1Left_04202017_Exp42/magnitude.vtk -scale 0 -type char -o  Kidney1Left_04202017_Exp42/roi.vtk
% $ c3d Kidney1Left_04202017_Exp42/magnitude.0001.vtk Kidney1Left_04202017_Exp42/roi.vtk -lstat
% LabelID        Mean        StdD         Max         Min       Count     Vol(mm^3)        Extent(Vox)
%     0           nan         nan  9368.31422     1.49276      327317    809045.156    256   256     5
%     1    6673.89057   752.07499  8121.13634  5553.91928         363       897.245     14    13     3
% $ c3d Kidney1Left_04202017_Exp42/magnitude.0001.vtk Kidney1Left_04202017_Exp42/magnitude.0003.vtk   -scale -1 -add Kidney1Left_04202017_Exp42/roi.vtk -lstat
% LabelID        Mean        StdD         Max         Min       Count     Vol(mm^3)        Extent(Vox)
%     0           nan         nan  1250.38738  -772.98665      327317    809045.156    256   256     5
%     1     -13.51417    43.71531   105.45898  -132.20173         363       897.245     14    13     3
% $ c3d Kidney1Left_04202017_Exp42/magnitude.0002.vtk Kidney1Left_04202017_Exp42/magnitude.0004.vtk   -scale -1 -add Kidney1Left_04202017_Exp42/roi.vtk -lstat
% LabelID        Mean        StdD         Max         Min       Count     Vol(mm^3)        Extent(Vox)
%     0           nan         nan  1122.43795  -959.95822      327317    809045.156    256   256     5
%     1     -11.02637    35.81174    87.49669  -144.63378         363       897.245     14    13     3
stdmag = 43.71 ; 

% https://www.tandfonline.com/doi/full/10.3109/02656736.2011.557028
snrfactor = stdmag/sqrt(2.)*sqrt(2.)/(2*pi*alpha   * gamma * field  *echotime  );

for jjj = 1 : length(outputPathList )
    outputPath = outputPathList{jjj}
    fullpathsplit =  strsplit(outputPath  ,'/');
    outputdir = fullpathsplit{ length(fullpathsplit)-1}
    mkdir(outputdir)

    fulldata = load(outputPath);
    %temperaturedata = fulldata.img.newdeltaT;
    temperaturedata = fulldata.newDeltaT;
    magnitudedata   = fulldata.datMag;
    r2data          = fulldata.datR2s;


    header.PixelSpacing   = fulldata.img.hdrEx.PixelSpacing;
    header.SliceThickness = fulldata.img.sliceThickness ;
    
    nsteps = size(temperaturedata ,4);
    
    for iii = 1:nsteps
        fileout = sprintf('%s/temperature.%04d.vtk',outputdir ,iii)
        saveVTKFile(fileout ,'temperature', temperaturedata(:,:,:,iii), header)
    end

    nsteps = size(magnitudedata,4);
    for iii = 1:nsteps
        fileout = sprintf('%s/magnitude.%04d.vtk',outputdir ,iii)
        saveVTKFile(fileout ,'magnitude', magnitudedata(:,:,:,iii), header)
    end

    for iii = 1:nsteps
        fileout = sprintf('%s/sigma.%04d.vtk',outputdir ,iii)
        saveVTKFile(fileout ,'sigma', snrfactor * (magnitudedata(:,:,:,iii)).^(-1), header)
    end
end





