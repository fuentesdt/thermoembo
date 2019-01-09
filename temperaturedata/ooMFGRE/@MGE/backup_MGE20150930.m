classdef MGE < handle
    
    properties (SetAccess=protected)
        %% Format specific properties
        path=''; % Path to exam
        numImages=0; % Number of images in exam directory
        imageOrder='Unknown'; % The order in which images come (GE DICOM only; Mag, Phase, Re, Im)
        format='Unknown' % Image format (DICOM (GE)/Bruker)
        imsPerIm=0; % Number of DICOM images associated with one true MR image (GE DICOM only)
        DICOMdic='/FUS4/data2/CJM/gitMATLAB/common/dict/gems-dicom-dictionary.txt' % DICOM dictionary
        imageIndices=0;
        
        %% Magnet specific scan parameters
        centerFrequency=0; % Magnet center frequency
        
        numSlices=0; % Number of slice in acqusition
        firstSliceLocation=0; % Location of first slice
        lastSliceLocation=0; %Location of last slice
        sliceThickness=0; % Slice thickness
        sliceSpacing=0; % Spacing between slice (center to center)
        distBetweenSlices=0; % Distance between slices (edge to edge)
        sliceOrientation=0; % Orientation of slices
        sliceCorner=0; %Corner of slices
        patientPosition=0; %Patient position (e.g. Head First Supine)
        
        hdrEx=0;

        acqWidth=0; % Width of acquisition matrix
        acqHeight=0; % Height of acquistion matrix
        displayWidth=0; % Width of display matrix
        displayHeight=0; % Height of display matrix
        
        acquisitionDuration=0 % Total time of acqusition
        numTimePoints=0; % Number of time points in the acquistion
        
        repetitionTime=0;
        fov=0;
        RBW=0;
        
        zip=0;
        
        flipAngle=0
        
        echoTrainLength=0; % Number of echoes
        echoTimes=0; % Echo times (for MFAMFGRE first echo train only)
        
        %% Data properties
        %cdata=0; % Complex data
        cdataRealInt=0;
        cdataImagInt=0;
        armaNumeratorList=0; % ARMA cumerator coefficients
        armaDenominatorList=0; % ARMA denominator coefficients
        armaRootsList=0; % Roots of ARMA denominator coefficients
        sortString='-amp'; % String defining how ARMA coefficients are ordered (default=descending by magnitude)
        processingROI=0; % ROI where ARMA coefficients are calculated
        validPixels=0; %arma coefficients that return valid values.
        
        %% Sequence specific parameters
        nFlip=1; % Number of flip angles
        flipDecrement=1
        nEchoShifts=1; % Number of echo shifts
        echoShiftFlag=0

    end
    
    properties (SetAccess=public)
        
        driftROI=0; % ROI use for drift correction (Must be withing processing ROI)
        driftPeak=1; % The peak to use when perfoming the drift correction
        
        model='stmcb' % ARMA solution model
        modelOrder=1; % Order of ARMA model
        gpu=0; % Use GPU or not
        
        threshold=100; % Singal threshold for performing ARMA processing
        numProcs=12; % Number of processors to use for paralelll CPU processing
        
        alpha=.01; % Temperature sensitivity coefficient (Fudged negative sign)
        deltaTPeak=1; % Peak to use when calculating change in temperature
        firstImage=3; % First (baseline) image to use for temperature imaging
        
        %% Background rmeoval algorithm properties
        BGRemOptions=struct('order',3,'nth',1e-6);

    end
    
    properties (Dependent)
        echoSpacing % Echo spacing
        sliceLocations % Vector of slice positions
        HBW %Half bandwidth
        roiIndx %Indices of processing roi
        r2star % r2star in 1/s
        t2star % t2star in ms
        ppm % resonance frequency in ppm
        amp % amplitude (complex)
        deltaT % deltaT (degC)
        %deltaT_MA % deltaT using moving average (Shifted echo Only)
        amp0 % ARMA calculated amplitude at TE=0
        minTEs % The minimum TE indexed by time point
        flipAngles % The flip angle indexed by time point
        modelSignal % Signal from reconstructed ARMA model.
        modelSignal2 % 
        ppmBGRem % PPM map with background removed
        memoryUse % Check RAM usage (kb)
        cdata % Complex Data (double)
        armaNumerator;
        armaDenominator;
        armaRoots;
%         armaDenominatorList;
%         armaRootsList;
        
    end
    
    methods
        %% Constructor Method
        function obj=MGE(MRscanPath)
            % CONSTRUCTOR/COPY METHOD: Creates/Copies MGE objects
            %  Accepts char, MGE, and struct as input
            %  MGE(char) char must be a path to raw Bruker or GE DICOM data
            %  MGE(struct) struct must contain parameters for generating signal. See /FUS4/data2/CJM/gitMATLAB/projects/mfgreScripts/inSilicoPhantoms/arma/genME.m
            %  MGE(MGE) MGE hard  copy the input object
            
            function [FTout]=fileType(fileList)
                switch 1
                    case ~isempty(regexp(fileList(1).name,('MR'), 'once'))
                        FTout='DICOM';
                    case ~isempty(regexp(fileList(1).name,('DCM'), 'once'))
                        FTout='DICOM';
                    otherwise
                        FTout='Bruker';
                end
            end %define function to determine filetype(GE/Bruker) if input is char.
            
            obj.path=MRscanPath;
            
            if ischar(MRscanPath)
                
                % Get a list of the files in the directory
                fileList=dir(MRscanPath);
                
                % Remove '.' and '..'
                fileList=fileList(3:numel(fileList));
                
                % Check/set format 
                obj.format=fileType(fileList);
                
                % Set properties/read data depending on format
                
                % For DICOM images
                if strcmp(obj.format,'DICOM');
                    
                    % Read first 4 headers,
                    hdr1=dicominfo([MRscanPath '/' fileList(1).name], 'dictionary',obj.DICOMdic);
                    hdr2=dicominfo([MRscanPath '/' fileList(2).name], 'dictionary',obj.DICOMdic);
                    hdr3=dicominfo([MRscanPath '/' fileList(3).name], 'dictionary',obj.DICOMdic);
                    hdr4=dicominfo([MRscanPath '/' fileList(4).name], 'dictionary',obj.DICOMdic);
                    
                    obj.hdrEx=hdr1;

                    % Determine number/type of images in folder (Mag, Phase,Real,Imag)
                    tmpdtype=unique([hdr1.RawDataType_ImageType,hdr2.RawDataType_ImageType,hdr3.RawDataType_ImageType,hdr4.RawDataType_ImageType]);
                    formatString={'Magnitude','Phase','Real','Imaginary'};
                    obj.imageOrder='';

                    % Set number of reconstructed images per actual image
                    obj.imsPerIm=numel(tmpdtype);
                    
                    % Set the order of reconstructed images
                    for nn=1:numel(tmpdtype)
                        obj.imageOrder=cat(2,formatString{tmpdtype+1});
                    end
                    
                    % Error check for insufficient data mus have real/imag or mag/phase
                    % TODO The code will take mag/phase but will break upon
                    % processing
                    %                     if numel(tmpdtype)<=1 || not(any(tmpdtype==1) && any(tmpdtype==0)) && not(any(tmpdtype==2) && any(tmpdtype==3))
                    %                         display('Wrong data types')
                    %                     end
                    
                    if tmpdtype==0;
                        display('Magnitude only; functionality limited')
                    elseif tmpdtype==[2 3];
                        display('Real + Imaginary')
                    end
                    
                    % Hardware Properties
                    obj.centerFrequency=hdr1.ImagingFrequency;
                    
                    obj.repetitionTime=hdr1.RepetitionTime;
                    obj.fov=hdr1.ReconstructionDiameter;
                    obj.RBW=hdr1.PixelBandwidth*128;
                    
                    % Rx Properties
                    obj.flipAngle=hdr1.FlipAngle;
                    obj.nFlip=hdr1.UserData11;
                    obj.flipDecrement=hdr1.UserData12;
                   
                    % Slice Rx Properties
                    obj.numSlices=double(hdr1.LocationsInAcquisition);
                    obj.sliceThickness=hdr1.SliceThickness;
                    obj.sliceSpacing=hdr1.SpacingBetweenSlices;
                    obj.firstSliceLocation=hdr1.FirstScanLocation;
                    obj.lastSliceLocation=hdr1.LastScanLocation;
                    obj.sliceOrientation=hdr1.ImageOrientationPatient;
                    obj.patientPosition=hdr1.PatientPosition;
                    
                    % Adjust slice thickness/locations if ZIP is on
                    if strcmp(hdr1.MRAcquisitionType,'3D')
                        obj.zip=round(obj.sliceThickness/obj.sliceSpacing);
                        obj.numSlices=(obj.zip*obj.numSlices);
                        obj.sliceThickness=hdr1.SliceThickness/obj.zip;
                        
                        if obj.lastSliceLocation>obj.firstSliceLocation && obj.zip>1
                            obj.lastSliceLocation=hdr1.LastScanLocation+obj.sliceThickness/2;
                             obj.firstSliceLocation=obj.firstSliceLocation-obj.sliceThickness/2;
                        elseif obj.lastSliceLocation<obj.firstSliceLocation
                            obj.lastSliceLocation=hdr1.LastScanLocation-obj.sliceThickness/2;
                            obj.firstSliceLocation=obj.firstSliceLocation+obj.sliceThickness/2;
                        end
                    end
                    
                    %Following block to account for ZIP2
%                     if strcmp(hdr1.MRAcquisitionType,'3D');
%                         obj.numSlices=(2*obj.numSlices);
%                         obj.sliceThickness=hdr1.SliceThickness/2;
%                         if obj.lastSliceLocation>obj.firstSliceLocation
%                             obj.lastSliceLocation=hdr1.LastScanLocation+obj.sliceThickness;
%                         elseif obj.lastSliceLocation<obj.firstSliceLocation
%                             obj.lastSliceLocation=hdr1.LastScanLocation-obj.sliceThickness;
%                         end
%                     end
                    
                    obj.displayWidth=double(hdr1.Width);
                    obj.displayHeight=double(hdr1.Height);
                    
                    obj.echoTrainLength=hdr1.NumberOfEchoes;
                    obj.nEchoShifts=hdr1.UserData18;
                    
                    obj.numImages=numel(fileList);
                    
                    obj.numTimePoints=double((2/obj.imsPerIm)*obj.numImages/obj.echoTrainLength/obj.numSlices/obj.imsPerIm);
                    
                    obj.acquisitionDuration=hdr1.AcquisitionDuration*1e-6;
 
                    % Preallocate temporary arrays
                    array1=zeros(obj.displayWidth,obj.displayHeight,obj.numTimePoints,obj.numSlices,obj.echoTrainLength,2);
                    time1=ones(obj.numSlices,obj.echoTrainLength,2);
                    tmpImageIndices=zeros(size(time1));
                    EchoTimes=zeros(1,obj.echoTrainLength);
                    tmpSliceLocations=zeros(obj.numSlices,3);
                    
                    % Loop over all images
                    for ii=1:obj.numImages
                        
                        % Print progress
                        if sum(ii==round([.05:.05:1]*obj.numImages));
                            disp([num2str(round(100*ii/obj.numImages)),'% Complete'])
                        end
                        
                        % Read header and find where each image belongs in the image stack
                        % TODO read phase directly header so you dont need to use the time1 array. 
                        % TODO avoid reading every header.
                        tmphdr=dicominfo([MRscanPath '/' fileList(ii).name], 'dictionary',obj.DICOMdic);
                        dataType=tmphdr.RawDataType_ImageType;
                        %SliceNumber=find(double(tmphdr.SliceLocation)==double(obj.sliceLocations));
                        SliceNumber=find((abs(double(tmphdr.SliceLocation)-double(obj.sliceLocations)))<10e-4);

                        % If data are real/imaginary; (e.g. rhrccrtl=28)
                        if ((dataType==2 || dataType==3) &&  ~isempty(SliceNumber))
                            
                            EchoTimes(ii)=tmphdr.EchoTime;
                            EchoNumber=tmphdr.EchoNumbers;
                            
                            array1(:,:,int16(time1(SliceNumber,EchoNumber,dataType-1)),SliceNumber,EchoNumber,dataType-1)=dicomread([MRscanPath '/' fileList(ii).name]);
                            time1(SliceNumber,EchoNumber,dataType-1)=time1(SliceNumber,EchoNumber,dataType-1)+1;
                            tmpSliceLocations(SliceNumber,:)=tmphdr.ImagePositionPatient;
                            tmpImageIndices(SliceNumber,EchoNumber,dataType-1)=ii;
                        
                        % If data are magnitude only; (e.g. rhrcctrl=17)
                        elseif (dataType==0 &&  ~isempty(SliceNumber))
                            
                            EchoTimes(ii)=tmphdr.EchoTime;
                            EchoNumber=tmphdr.EchoNumbers;
                            
                            array1(:,:,int16(time1(SliceNumber,EchoNumber,dataType+1)),SliceNumber,EchoNumber,1)=dicomread([MRscanPath '/' fileList(ii).name]);
                            array1(:,:,int16(time1(SliceNumber,EchoNumber,dataType+1)),SliceNumber,EchoNumber,2)=zeros(obj.displayWidth,obj.displayHeight);
                            time1(SliceNumber,EchoNumber,dataType+1)=time1(SliceNumber,EchoNumber,dataType+1)+1;
                            tmpSliceLocations(SliceNumber,:)=tmphdr.ImagePositionPatient;
                            tmpImageIndices(SliceNumber,EchoNumber,dataType+1)=ii;
                        end
                        
                    end

                    obj.imageIndices=tmpImageIndices;
                    obj.sliceCorner=tmpSliceLocations;

                    % Permute image stack into row,col,slice,echo,timepoint format
                    obj.cdataRealInt=permute(int16(array1(:,:,:,:,:,1)),[1 2 4 5 3]);
                    obj.cdataImagInt=permute(int16(array1(:,:,:,:,:,2)),[1 2 4 5 3]);
                    
                    % Set echo time values
                    obj.echoTimes=nonzeros(unique(EchoTimes));     
                    
                    % If the acqusition is in 3D apply fftshift in slice direction to fix GE bug.
                    if strcmp(hdr1.MRAcquisitionType,'3D') && ~strcmp(hdr1.SeriesDescription,'NOT DIAGNOSTIC: Custom Recon')
 
                        tempdatamat=(ifft(fftshift(fft(complex(double(obj.cdataRealInt),double(obj.cdataImagInt)),[],3),3),[],3));
                        obj.cdataRealInt=real(tempdatamat);
                        obj.cdataImagInt=imag(tempdatamat);
                    end
                    
                % For Bruker Raw data    
                elseif strcmp(obj.format,'Bruker');
                    
                    % Read image data
                    [~,img,~,~] = Bruker_MGE(MRscanPath);
                    
                    % Read header data
                    hdr = readBrukerHeader([MRscanPath 'method']);
                    hdr2 = readBrukerHeader([MRscanPath '/pdata/1/visu_pars']);
                    
                    obj.cdataRealInt=int16(real(img));
                    obj.cdataImagInt=int16(imag(img));
                    
                    obj.echoTimes=hdr.EffectiveTE;
                    obj.echoTrainLength=numel(img(1,1,1,:,1));
                    
                    obj.centerFrequency=hdr2.VisuAcqImagingFrequency;
                    
                    obj.displayWidth=numel(img(1,:,1,1,1));
                    obj.displayHeight=numel(img(:,1,1,1,1));
                    
                    % TODO Check the geometry for 2D Scans in slice direction
                    obj.numSlices=hdr.PVM_SPackArrNSlices;
                    obj.numSlices=numel(img(1,1,:,1,1));
                    obj.sliceThickness=hdr.PVM_SliceThick;
                    obj.sliceSpacing=hdr.PVM_SPackArrSliceGap;
                    obj.firstSliceLocation=0;min(hdr.PVM_FovCm);
                    obj.lastSliceLocation=(obj.numSlices-1)*obj.sliceSpacing;
                    obj.distBetweenSlices=hdr.PVM_SPackArrSliceDistance;
                    
                    obj.numImages=numel(img(1,1,1,1,:))*numel(img(1,1,1,:,1))*numel(img(1,1,:,1,1));
                    obj.imageOrder='NA';
                    obj.imsPerIm=1;
                    
                    obj.acqWidth=hdr.PVM_EncMatrix(2);
                    obj.acqHeight=hdr.PVM_EncMatrix(1);       
                    
                    obj.flipAngle=hdr.PVM_ExcPulseAngle;
                    
                    obj.numTimePoints=obj.numImages/obj.echoTrainLength/obj.numSlices/obj.imsPerIm;
                    % TODO Find acqusition time in Bruker header
                    tmpacqDur=hdr.PVM_ScanTimeStr; %example 0h6m18s0ms
                    hourindx=strfind(tmpacqDur,'h');
                    minindx=strfind(tmpacqDur,'m');
                    secindx=strfind(tmpacqDur,'s');
                    msecindx=strfind(tmpacqDur,'ms');
                    obj.acquisitionDuration=3600*str2double(tmpacqDur(1:(hourindx-1)))+60*str2double(tmpacqDur((hourindx+1):(minindx-1)))+str2double(tmpacqDur((minindx+1):(secindx-1)))+(10^-3)*str2double(tmpacqDur((secindx+1):(msecindx-1)));
                end
                
            % Hardcopy input if class(input)=MGE    
            elseif isa(MRscanPath,'MGE')
                propsCell=properties('MGE');
                numProps=numel(propsCell);
                for ii=1:numProps
                    isDep=findprop(MRscanPath,propsCell{ii});
                    if isDep.Dependent==0
                        eval(['obj.' propsCell{ii} '=MRscanPath.' propsCell{ii},';']);
                    end
                end

                % Create simulation data from structure
            elseif isstruct(MRscanPath)
                [tmpOut] = genME( MRscanPath );
                obj.path='simulation';
                
                obj.cdataRealInt=real(tmpOut.cdata);
                obj.cdataImagInt=imag(tmpOut.cdata);
                %tmpOut.tdata=permute(wonoise,[1 2 4 3]);
                obj.echoTrainLength=numel(tmpOut.EchoTimes);
                obj.echoTimes=tmpOut.EchoTimes;
                obj.centerFrequency=tmpOut.hdr.ImagingFrequency;
                obj.displayWidth=tmpOut.hdr.Width;
                obj.displayHeight=tmpOut.hdr.Height;
                obj.numSlices=tmpOut.hdr.LocationsInAcquisition;
                obj.numTimePoints=1;
               
            else
                disp('Incompatible input')
            end
            
            
            obj.cdataRealInt=reshape(permute(obj.cdataRealInt,[1 2 3 5 4]),[],obj.echoTrainLength);
            obj.cdataImagInt=reshape(permute(obj.cdataImagInt,[1 2 3 5 4]),[],obj.echoTrainLength);

        end
        
        %% Set Methods
        
        function setProcessingROI(obj,varargin)
            % SET PROCESSING ROI: Used to define the ROI where ARMA
            % coefficients will be calculated. 
            % 
            % No argument: Displays mddle slice for interactive ROI
            % placement
            %
            % One argument: 4 element vector defining rectangular ROI coordinates in (x,y,w,h) form
            if nargin==1
                disp('Select ROI for processing')
                imagesc(abs(obj.cdata(:,:,round(obj.numSlices/2),1,1)))
                obj.processingROI=round(getrect);
                close
            elseif numel(varargin{1})==1
                disp('Select center of ROI for processing')
                imagesc(abs(obj.cdata(:,:,round(obj.numSlices/2),1,1)))
                [xx yy]=getpts;
                xx=round(xx);
                yy=round(yy);
                obj.processingROI=[xx-round(varargin{1}/2) yy-round(varargin{1}) varargin{1} varargin{1}];
                close
            elseif numel(varargin{1})==4
                obj.processingROI=varargin{1};
            end
            
            obj.calcArmaCoeffs
        end
        
        function setDriftROI(obj,varargin)
            % SET DRIFT CORRECTION ROI: Used to define the ROI that will be
            % used to correct for B0 field drift. The drift ROI must be
            % completely contained by the processing ROI.
            %
            % No argument: Displays mddle slice for interactive ROI
            % placement
            %
            % One argument: 4 element vector defining rectangular ROI coordinates in (x,y,w,h) form
            if nargin==1
                disp('Select ROI for processing')
                imagesc(abs(obj.cdata(:,:,round(obj.numSlices/2),1,1)))
                obj.driftROI=round(getrect);
                close
            else
                obj.driftROI=varargin{1};
            end
        end
        
        function setSortString(obj,newSortFlag)
            % SET SORTING STRING: Sorts the ARMA peaks based on
            % string input
            % 
            % Accepted inputs:
            % '-(+)amp':peaks sorted descending(ascending) by amplitude
            % '-(+)r2star':peaks sorted descending(ascending) by r2*
            % '-(+)t2star':peaks sorted descending(ascending) by t2*
            % '-(+)ppm':peaks sorted descending(ascending) by ppm
            
            function sortCoeffs(obj)
                
                switch 1
                    case strcmp(obj.sortString,'amp')||strcmp(obj.sortString,'+amp')
                        
                        pp=zeros(size(obj.armaRootsList)); qq=pp;
                        for ii=1:obj.modelOrder
                            pp=pp+repmat(obj.armaNumeratorList(:,ii),[1 obj.modelOrder]).*(obj.armaRootsList).^(obj.modelOrder-ii);
                            qq=qq+(obj.modelOrder+1-ii)*repmat(obj.armaDenominatorList(:,ii),[1 obj.modelOrder]).*(obj.armaRootsList).^(obj.modelOrder-ii);
                        end
                        tmpamp=pp./qq;
                        
                        sortarray=tmpamp;
                        sortstring='ascend';
                    case strcmp(obj.sortString,'-amp')
                        
                        pp=zeros(size(obj.armaRootsList)); qq=pp;
                        for ii=1:obj.modelOrder
                            pp=pp+repmat(obj.armaNumeratorList(:,ii),[1 obj.modelOrder]).*(obj.armaRootsList).^(obj.modelOrder-ii);
                            qq=qq+(obj.modelOrder+1-ii)*repmat(obj.armaDenominatorList(:,ii),[1 obj.modelOrder]).*(obj.armaRootsList).^(obj.modelOrder-ii);
                        end
                        tmpamp=pp./qq;
                        
                        sortarray=tmpamp;
                        sortstring='descend';
                    case strcmp(obj.sortString,'ppm')||strcmp(obj.sortString,'+ppm')
                        sortarray=-imag(log(obj.armaRootsList))/2/pi/(obj.echoSpacing*10^-3*obj.centerFrequency);
                        sortstring='ascend';
                    case strcmp(obj.sortString,'-ppm')
                        sortstring='descend';
                        sortarray=-imag(log(obj.armaRootsList))/2/pi/(obj.echoSpacing*10^-3*obj.centerFrequency);
                    case strcmp(obj.sortString,'r2s')||strcmp(obj.sortString,'+r2s')
                        sortstring='descend';
                        sortarray=-1000*real(log((obj.armaRootsList))./obj.echoSpacing);
                    case strcmp(obj.sortString,'-r2s')
                        sortarray=-1000*real(log((obj.armaRootsList))./obj.echoSpacing);
                        sortstring='ascend';
                    case strcmp(obj.sortString,'t2s')||strcmp(obj.sortString,'+t2s')
                        sortarray=1/(-1000*real(log((obj.armaRootsList))./obj.echoSpacing));
                        sortstring='ascend';
                    case strcmp(obj.sortString,'-t2s')
                        sortarray=1/(-1000*real(log((obj.armaRootsList))./obj.echoSpacing));
                        sortstring='descend';
                end
                
                [~,sortindx]=sortlind(nan2zero(sortarray),2,sortstring);
                
                obj.armaNumeratorList=obj.armaNumeratorList(sortindx);
                tmpArray=obj.armaDenominatorList(:,2:end);
                tmpArray=tmpArray(sortindx);
                obj.armaDenominatorList(:,2:end)=tmpArray;
                obj.armaRootsList=obj.armaRootsList(sortindx);
                
            end
            
            obj.sortString=newSortFlag;
            sortCoeffs(obj)
        end
        
        function cropSlices(obj,sliceVec)
            tmpdata=obj.cdata(:,:,sliceVec,:,:);
            tmpsliceLocs=obj.sliceLocations(sliceVec);
            obj.firstSliceLocation=tmpsliceLocs(1);
            obj.lastSliceLocation=tmpsliceLocs(end);
            obj.numSlices=numel(sliceVec);
            obj.setCdata(tmpdata);
            
        end
        
        function cropTimePoints(obj,TPVec)
            tmpdata=obj.cdata(:,:,:,:,TPVec);
            obj.numTimePoints=numel(TPVec); 
            obj.setCdata(tmpdata);
            %obj.cdataRealInt=int16(obj.cdataRealInt(:,:,:,:,TPVec));
            %obj.cdataImagInt=int16(obj.cdataImagInt(:,:,:,:,TPVec));
            %obj.calcArmaCoeffs;
        end
        
        function cropEchoes(obj,echoVec)
            tmpdata=obj.cdata(:,:,:,echoVec,:);
            obj.echoTrainLength=numel(echoVec);
            obj.echoTimes=obj.echoTimes(echoVec);
            %obj.cdata=obj.cdata(:,:,:,echoVec,:);
            obj.setCdata(tmpdata)
%           obj.cdataRealInt=obj.cdataRealInt(:,:,:,echoVec,:);
%           obj.cdataImagInt=obj.cdataImagInt(:,:,:,echoVec,:);
            %obj.calcArmaCoeffs
        end     
        
        function setCdata(obj,newCdata)
            obj.cdataRealInt=reshape(permute(real(newCdata),[1 2 3 5 4]),[],obj.echoTrainLength);
            obj.cdataImagInt=reshape(permute(imag(newCdata),[1 2 3 5 4]),[],obj.echoTrainLength);
        end
        
        %% Dependent Property Methods
        function cdata=get.cdata(obj)
            
            cdata=double(permute(reshape(complex(obj.cdataRealInt,obj.cdataImagInt),obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.echoTrainLength),[1 2 3 5 4]));
            
        end
        
        function armaNumerator=get.armaNumerator(obj)
            
            tmpArray=nan((obj.processingROI(4)+1)*(obj.processingROI(3)+1)*double(obj.numSlices)*double(obj.numTimePoints),obj.modelOrder);
            tmpArray(obj.validPixels,:)=obj.armaNumeratorList;
            
            tmpArray2=reshape(tmpArray,[obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder]);
            
            armaNumerator=zeros(obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.modelOrder);
            armaNumerator(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=tmpArray2;
        end
        
        function armaDenominator=get.armaDenominator(obj)
            %armaDenominator=reshape(obj.armaDenominatorList,[obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder+1]);
            
            tmpArray=nan((obj.processingROI(4)+1)*(obj.processingROI(3)+1)*double(obj.numSlices)*double(obj.numTimePoints),obj.modelOrder+1);
            tmpArray(obj.validPixels,:)=obj.armaDenominatorList;
            
            tmpArray2=reshape(tmpArray,[obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder+1]);
            
            armaDenominator=zeros(obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.modelOrder+1);
            armaDenominator(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=tmpArray2;
            
        end
        
        function armaRoots=get.armaRoots(obj)
            
            tmpArray=nan((obj.processingROI(4)+1)*(obj.processingROI(3)+1)*double(obj.numSlices)*double(obj.numTimePoints),obj.modelOrder);
            tmpArray(obj.validPixels,:)=obj.armaRootsList;
            
            tmpArray2=reshape(tmpArray,[obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder]);
            
            armaRoots=zeros(obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.modelOrder);
            armaRoots(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=tmpArray2;
            
        end
        
        function r2star=get.r2star(obj)
            tmpr2star=-1000*real(log((obj.armaRootsList))./obj.echoSpacing);
            
            tmpr22=nan((obj.processingROI(4)+1)*(obj.processingROI(3)+1)*obj.numSlices.*obj.numTimePoints,obj.modelOrder);
            tmpr22(obj.validPixels,:)=tmpr2star;
            
            tmpr23=reshape(tmpr22,[obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder]);
            
            r2star=zeros(obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.modelOrder);
            r2star(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=inf2zero(tmpr23);
        end
        
        function t2star=get.t2star(obj)
            tmpt2star=-obj.echoSpacing./real(log(obj.armaRoots));
            t2star=zeros(obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.modelOrder);
            t2star(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=tmpt2star;
        end
        
        function ppm=get.ppm(obj)
            
            tmpppm=-imag(log(obj.armaRootsList))/2/pi/(obj.echoSpacing*10^-3*obj.centerFrequency);
            
            tmpppm2=nan((obj.processingROI(4)+1)*(obj.processingROI(3)+1)*obj.numSlices.*obj.numTimePoints,obj.modelOrder);
            tmpppm2(obj.validPixels,:)=tmpppm;
            
            tmpppm3=reshape(tmpppm2,[obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder]);
            
            ppm=zeros(obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.modelOrder);
            ppm(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=tmpppm3;
        end
        
        function roiIndx=get.roiIndx(obj)
            roiIndx=ROI2indx(obj.processingROI);
        end
        
        function echoSpacing=get.echoSpacing(obj)
            echoSpacing=obj.echoTimes(2)-obj.echoTimes(1);
            if numel(unique(echoSpacing))~=1
                disp('Echoes are not evenly spaced')
            end
        end
        
        function sliceLocations=get.sliceLocations(obj)
            sliceLocations=linspace(obj.firstSliceLocation,obj.lastSliceLocation,obj.numSlices);
        end
        
        function HBW=get.HBW(obj)
            HBW=1/(2*obj.echoSpacing*obj.centerFrequency*10^-3);
        end
        
        function ppmBGRem=get.ppmBGRem(obj)
            tmpppm=imcrop2(obj.ppm,obj.processingROI);
            tmpamp=nan2zero(imcrop2(abs(obj.amp).^2,obj.processingROI));

            BGrem=(obj.HBW/pi)*bgsup3((pi/obj.HBW)*tmpppm,tmpamp,obj.BGRemOptions);
            
            ppmBGRem=zeros(obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.modelOrder);
            ppmBGRem(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=BGrem;
        end
        
        function amp=get.amp(obj)
            
            pp=zeros(size(obj.armaRootsList)); qq=pp;
            for ii=1:obj.modelOrder
                pp=pp+repmat(obj.armaNumeratorList(:,ii),[1 obj.modelOrder]).*(obj.armaRootsList).^(obj.modelOrder-ii);
                qq=qq+(obj.modelOrder+1-ii)*repmat(obj.armaDenominatorList(:,ii),[1 obj.modelOrder]).*(obj.armaRootsList).^(obj.modelOrder-ii);
            end
            tmpamp=pp./qq;
            
            tmpamp2=nan((obj.processingROI(4)+1)*(obj.processingROI(3)+1)*obj.numSlices.*obj.numTimePoints,obj.modelOrder);
            tmpamp2(obj.validPixels,:)=tmpamp;
            
            tmpamp3=reshape(tmpamp2,[obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder]);
            
            amp=zeros(obj.displayHeight,obj.displayWidth,obj.numSlices,obj.numTimePoints,obj.modelOrder);
            amp(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=tmpamp3;
        end
        
        function amp0=get.amp0(obj)
            minTEmat=repmat(permute(obj.minTEs,[1 4 3 2]),[numel(obj.cdata(:,1,1,1,1)) numel(obj.cdata(1,:,1,1,1)) numel(obj.cdata(1,1,:,1,1)) 1, obj.modelOrder]);
            amp0=obj.amp.*exp(2.* pi.* 1i .* obj.ppm .* minTEmat .* obj.centerFrequency .* 10^-3 ) .* exp((minTEmat).*obj.r2star./1000);
        end
        
        function minTEs=get.minTEs(obj)
%             tmpminTEs=linspace(obj.echoTimes(1),obj.echoTimes(1)+obj.echoSpacing,obj.nEchoShifts+1);
%             tmpminTEs=tmpminTEs(1:end-1);
%             minTEs=repmat(tmpminTEs,[1,ceil(obj.numTimePoints/obj.nEchoShifts)]);
%             minTEs=minTEs(1:obj.numTimePoints);
        minTEs=obj.echoTimes(1)*ones(1,obj.numTimePoints);

        end
        
        function flipAngles=get.flipAngles(obj)
            % This is convoluted, fix for the sake of
            % simplicity/readability
            tmpAlphas=repmat(obj.flipAngle:-obj.flipDecrement:obj.flipAngle-(obj.nFlip-1)*obj.flipDecrement,[1,obj.numTimePoints/obj.nEchoShifts/obj.nFlip]);
            tmpAlphas2=ones(1,obj.numTimePoints);
            for ii=1:obj.nEchoShifts
                tmpAlphas2(ii:obj.nEchoShifts:end)=tmpAlphas;
            end
            flipAngles=tmpAlphas2;
        end
        
        function deltaT=get.deltaT(obj)
            
            %uncorrDeltaT=-(1/obj.alpha)*(obj.ppm(:,:,:,:,obj.deltaTPeak)-repmat(obj.ppm(:,:,:,obj.firstImage,obj.deltaTPeak),[1 1 1 numel(obj.ppm(1,1,1,:,1)) 1]));
            
            uwppm=(obj.HBW/pi)*(unwrap((pi/obj.HBW)*(obj.ppm(:,:,:,:,obj.deltaTPeak)),[],4));
            baseline=repmat(uwppm(:,:,:,obj.firstImage,1),[1 1 1 numel(obj.ppm(1,1,1,:,1)) 1]);
            uncorrDeltaT=-(1/obj.alpha)*(uwppm-baseline);
            %uncorrDeltaT=-(1/obj.alpha)*(uwppm-repmat(uwppm(:,:,:,obj.firstImage,obj.deltaTPeak),[1 1 1 numel(obj.ppm(1,1,1,:,1)) 1])),[],4));
            
            %uncorrDeltaT=-(1/obj.alpha)*(obj.HBW/pi)*(unwrap((pi/obj.HBW)*(obj.ppm(:,:,:,:,obj.deltaTPeak)-repmat(obj.ppm(:,:,:,obj.firstImage,obj.deltaTPeak),[1 1 1 numel(obj.ppm(1,1,1,:,1)) 1])),[],4));
            
            
            switch 1
                case obj.numTimePoints<2
                    disp('Data does not have multiple timepoints')
                    deltaT=zeros(size(obj.ppm));
                case isequal(obj.driftROI,0) || isequal(obj.driftROI,[0 0 0 0]);
                    % obj.deltaT=(1/obj.alpha)*obj.ppm(:,:,:,:,obj.deltaTPeak)-obj.ppm(:,:,:,obj.firstImage,obj.deltaTPeak);
                    deltaT=uncorrDeltaT;
                case ~isequal(obj.driftROI,0)
                    resizeVec=size(obj.ppm);
                    resizeVec(4)=1;
                    tmpindx=ROI2indx(obj.driftROI);
                    driftCorrect=repmat(mean(mean(uncorrDeltaT(tmpindx(1):tmpindx(2),tmpindx(3):tmpindx(4),round(obj.numSlices/2),:),2),1),resizeVec);
                    %driftCorrect=repmat(mean(mean(imcrop2((1/obj.alpha)*obj.ppm(:,:,:,:,obj.driftPeak)-repmat(obj.ppm(:,:,:,obj.firstImage,obj.driftPeak),[1 1 1 ])),2),1),resizeVec);
                    deltaT=uncorrDeltaT-driftCorrect(:,:,:,:,1);
            end
            
        end
        
        function modelSignal=get.modelSignal(obj)
            TEmat=zeros(obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.echoTrainLength,obj.numTimePoints,obj.modelOrder);
            TEmat(:,:,:,1,:,:)=repmat(permute(obj.minTEs,[1 5 3 4 2 6]),[obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,1,1,obj.modelOrder]);
            
            for ii=2:obj.echoTrainLength
                TEmat(:,:,:,ii,:)=TEmat(:,:,:,1,:)+(double(ii)-1)*obj.echoSpacing;
            end
            
            addEchoDim=@(pp) repmat(permute(pp(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:),[1 2 3 7 4 5 6]),[1 1 1 obj.echoTrainLength 1]);
            
            tmpModelSignal =  addEchoDim(obj.amp0) .* exp(2.* pi.* -1i .* addEchoDim(obj.ppm) .* TEmat .* obj.centerFrequency .* 10^-3 ) .* exp( -TEmat ./ addEchoDim(obj.t2star) );
            
            modelSignal=zeros([size(obj.cdata),obj.modelOrder]);
            modelSignal(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:,:)=tmpModelSignal;
                       
        end
        
        function modelSignal2=get.modelSignal2(obj)
            N=obj.echoTrainLength;
            modelSignal2=zeros(size(obj.cdata));
            filtfun=@(xx,yy,nn) filter(xx,yy,[1 zeros(1,nn-1)]);
            tmpModelSig=zeros(numel(obj.armaNumeratorList(:,1)),N);
            
            for ii=1:numel(obj.armaNumeratorList(:,1))
                if abs(obj.armaDenominatorList(ii,:))>0;
                    tmpModelSig(ii,:)=filtfun(obj.armaNumeratorList(ii,:),obj.armaDenominatorList(ii,:),N);
                end
            end
            
            tmparray=reshape(tmpModelSig,obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.echoTrainLength,obj.numTimePoints);
            modelSignal2(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:)=tmparray;
        end
        
        function memoryUse=get.memoryUse(obj)
            memoryUse=0;
            propsCell=properties('MGE');
            numProps=numel(propsCell);
            for ii=1:numProps
                isDep=findprop(obj,propsCell{ii});
                if isDep.Dependent==0
                    tmpProp=eval(['obj.' propsCell{ii}]);
                    tmpProp2=whos('tmpProp');
                    tmpProp3=tmpProp2.bytes;
                    memoryUse=memoryUse+tmpProp3;
                end
            end
            memoryUse=memoryUse*1e-6
        end
        
        %% Misc Methods
        
        function calcArmaCoeffs(obj)
            % TODO Fix waitbars/add timers
            % TODO Make thresholding more efficient using mask/indexing
            
            clear obj.armaNumerator obj.armaDenominator obj.armaRoots
            
            % Reshape complex data
            if obj.processingROI==0
                obj.processingROI=[1 1 obj.displayHeight-1 obj.displayWidth-1];
            end
            
            tmpcdata=reshape(permute(obj.cdata(obj.roiIndx(1):obj.roiIndx(2),obj.roiIndx(3):obj.roiIndx(4),:,:,:),[1 2 3 5 4]),[],obj.echoTrainLength);
            
            procPixIdx=find(abs(tmpcdata(:,1))>obj.threshold & abs(real(tmpcdata(:,1)))~=32767 & abs(imag(tmpcdata(:,1)))~=32767);
            %procPixIdx=(1:numel(tmpcdata(:,1)))';
            procCdata=tmpcdata(procPixIdx,:);
            
            procNum=zeros(numel(procCdata(:,1)),obj.modelOrder);
            procDen=zeros(numel(procCdata(:,1)),obj.modelOrder+1);
            procRoots=procNum;
            
            switch 1
                case obj.numProcs==1;
                    for ii=1:numel(procCdata(:,1))
                        [procNum(ii,:),procDen(ii,:),procRoots(ii,:)]=csi_armaMGE(squeeze(procCdata(ii,:)),obj.modelOrder,obj.model);
                    end
                    
                case obj.numProcs>1;
%                     if matlabpool('size')==0
%                         % Fic this line to open pool of size numProcs
%                         % matlabpool obj.numProcs
%                     end

                    tmpModelOrder=obj.modelOrder;
                    tmpModel=obj.model;
                    
                    parfor (ii=1:numel(procCdata(:,1)),obj.numProcs)
                        tmpcdata2=squeeze(procCdata(ii,:));
                        
                        [procNum(ii,:),procDen(ii,:),procRoots(ii,:)]=csi_armaMGE(tmpcdata2,tmpModelOrder,tmpModel);
                        
                    end
                    
                    % TODO Incorporate GPU processing. Overload to keep data in GPU arrays 
                    
            end
           
            tmpNum=zeros(numel(tmpcdata(:,1)),obj.modelOrder);
            tmpDen=zeros(numel(tmpcdata(:,1)),obj.modelOrder+1);
            tmpRoots=tmpNum;
            
            tmpNum(procPixIdx,:)=procNum;
            tmpDen(procPixIdx,:)=procDen;
            tmpRoots(procPixIdx,:)=procRoots;
      
            validPix=find((abs(tmpDen(:,1))==1).*(abs(tmpNum(:,1))>1));
            
            obj.validPixels=validPix;
            obj.armaNumeratorList=tmpNum(validPix,:);
            obj.armaDenominatorList=tmpDen(validPix,:);
            obj.armaRootsList=tmpRoots(validPix,:);
            
        end
        % TODO Store coefficients as single precision?
 
        function fixFirstEcho(obj)
            close all
            %imagesc(abs(obj.cdata(:,:,1,1,1)))
            %[rr,cc]=getpts;
            rr=160; %round(rr);
            cc=113; %round(cc);
            plot(unwrap(angle(squeeze(obj.cdata(rr,cc,1,:,1)))),'k'); hold on;
            endflag=0;
            tmpcdata=obj.cdata;
            
            while endflag==0;
            nn=input('N?');
            if isnumeric(nn)
                if exist('htmp')
                delete(htmp)
                end
                tmpcdata(:,:,:,1,:)=obj.cdata(:,:,:,1,:).*exp(-2*pi*1i*2*obj.HBW*nn);
                htmp=plot(unwrap(angle(squeeze(tmpcdata(rr,cc,1,:,1)))),'r');
            elseif ischar(nn)
                endflag=1;
            end
            
            obj.setCdata(tmpcdata)
            end
                
        end
        
        function obj2=MGEechoShift(obj)
            
            obj2=MGE(obj);
            
            obj2.echoTrainLength=obj.echoTrainLength*obj.nEchoShifts;
            
            obj2.echoTimes=linspace(obj.echoTimes(1),obj.echoTimes(end)+(obj.echoSpacing-(obj.echoSpacing/obj.nEchoShifts)),obj2.echoTrainLength);
            
            tmpTimePoints=double(obj.numTimePoints)/double(obj.nEchoShifts);
            if double(int16(tmpTimePoints))==tmpTimePoints
                obj2.cropTimePoints([1:tmpTimePoints]);
            else
                obj2.cropTimePoints([1:floor(tmpTimePoints)])
            end
            
            obj2.nEchoShifts=1;
            
%             tmpSize=size(obj2.cdata);
%             tmpSize(4)=obj.nEchoShifts*obj.echoTrainLength;
%             tmpCdataRealInt=zeros(tmpSize,'int16');
%             tmpCdataImagInt=zeros(tmpSize,'int16');
            
            tmpcdata=zeros(obj.displayWidth,obj.displayHeight,obj.numSlices,obj2.echoTrainLength,obj2.numTimePoints);

            kk=1;
            for ii=1:obj2.numTimePoints;
                for jj=1:obj.nEchoShifts;
%                     tmpCdataRealInt(:,:,:,jj:obj.nEchoShifts:end,ii)=obj.cdataRealInt(:,:,:,:,kk);
%                     tmpCdataImagInt(:,:,:,jj:obj.nEchoShifts:end,ii)=obj.cdataImagInt(:,:,:,:,kk);
                    tmpcdata(:,:,:,jj:obj.nEchoShifts:end,ii)=obj.cdata(:,:,:,:,kk);
                    kk=kk+1;
                end
            end
%             obj2.cdataRealInt=tmpCdataRealInt;
%             obj2.cdataImagInt=tmpCdataImagInt;
            obj.setCdata(tmpcdata)
            
            obj2.echoShiftFlag=1;
            obj2.calcArmaCoeffs;
        end
        
        function obj2=applyMAfilter(obj)
            
            obj2=MGE(obj);
            obj2.cdataRealInt=zeros(obj2.displayWidth,obj2.displayHeight,obj2.numSlices,obj2.echoTrainLength*obj2.nEchoShifts,obj2.numTimePoints,'int16');
            obj2.cdataImagInt=zeros(obj2.displayWidth,obj2.displayHeight,obj2.numSlices,obj2.echoTrainLength*obj2.nEchoShifts,obj2.numTimePoints,'int16');
            obj2.armaNumerator=zeros(obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder);
            obj2.armaDenominator=zeros(obj.processingROI(4)+1,obj.processingROI(3)+1,obj.numSlices,obj.numTimePoints,obj.modelOrder+1);
            obj2.armaRoots=obj.armaNumerator;
            
            %cropVec=[zeros(1,obj.nEchoShifts*obj.numTimePoints)];
            tmpCropVec=[1:obj.nEchoShifts];
            cropVec=zeros(size(tmpCropVec),'int16');
            
            for ii=1:obj.numTimePoints
                %startInd=[1:obj.nEchoShifts:obj.nEchoShifts*obj.numTimePoints];
                if lt(ii,obj.nEchoShifts)
                    cropVec(1:obj.nEchoShifts)=[1 2 3 4];
                else
                    cropVec(1:obj.nEchoShifts)=[tmpCropVec];
                    tmpindx=mod(ii,obj.nEchoShifts)+1;
                    tmpCropVec(tmpindx)=tmpCropVec(tmpindx)+obj.nEchoShifts;
                end
                tmpObj=MGE(obj);
                tmpObj.cropTimePoints(cropVec)
                tmpObj2=tmpObj.MGEechoShift;
                obj2.cdataRealInt(:,:,:,:,ii)=tmpObj2.cdataRealInt;
                obj2.cdataImagInt(:,:,:,:,ii)=tmpObj2.cdataImagInt;
                obj2.armaNumerator(:,:,:,ii,:)=tmpObj2.armaNumerator;
                obj2.armaDenominator(:,:,:,ii,:)=tmpObj2.armaDenominator;
                obj2.armaRoots(:,:,:,ii,:)=tmpObj2.armaRoots;
                obj2.echoTrainLength=obj.echoTrainLength*obj.nEchoShifts;
                obj2.nEchoShifts=1;
                obj2.echoTimes=tmpObj2.echoTimes;
            end
            obj2.firstImage=obj.nEchoShifts;
            
        end
        
    end
    
end