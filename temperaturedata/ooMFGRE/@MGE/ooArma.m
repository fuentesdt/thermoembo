function [data,roi]=ooArma(images,varargin)
%% Performs ARMA processing to give spectral information.

% Usage: [ims2]=arma(ims,'OptionalParamString','OptionalParamValue'....)]

% Outputs
% 'ims2': structure containing spectral information
%     'ims2.t2star': t2star information [ms] 
%     'ims2.ppm': ppm information [ppm]
%     'ims2.amp': amplitude information [AU, complex]
%     'ims2.t1': t1 information (work in progress) [ms]
%     'ims2.ROI': ROI for processing [4 element vector]
%     'ims2.droi': ROI used for drift calculation [4 element vector]
%     'ims2.drift': temperature correction from drift [degC]
%     'ims2.deltaT': change in temperature [degC]
% 'roi': ROI used for processing

% Required inputs:
% 'ims': Imaging data. This must be the structure output of importMFGRE.m

% Optional inputs:
% 'numPeaks': Number of peaks desired [positive integer]
% 'alpha': Temperature sensitivity coefficient [ppm/degC]
% 'firstImage': First image used for temperature calculation [positive integer]
% 'ROI': ROI for processing [4 element vector]
% 'driftROI': ROI used for drift calculation [4 element vector]
% 'thresh': Minimum threshold needed to perform processing
% 'procs': Number of CPUs to use
% 'model': Type of model used analyze signal
%   options: 'stmcb'(Steiglitz-McBride,default),'prony'(Prony)
% 'gpu': Use GPU (*Prony only)
%   options: '0'(GPU off,default),'1'(GPU on)
% 'sort': String to indicate how peaks should be sorted
%   options: 'amp' (amplitude ascending) 'r2s' (r2s ascending) etc.


%% Parse inputs using input parser class
inPars=inputParser;
addRequired(inPars,'images',@isstruct);
addOptional(inPars,'numPeaks',1,@isnumeric)
addOptional(inPars,'alpha',-.01,@isnumeric)
addOptional(inPars,'firstImage',3,@isnumeric)
addOptional(inPars,'ROI', 0,@isnumeric)
addOptional(inPars,'driftROI', 0,@isnumeric)
addOptional(inPars,'thresh', 100,@isnumeric)
addOptional(inPars,'procs', 0,@isnumeric)
% TODO can we clean up the flags to use function handles ? 
addOptional(inPars,'model','stmcb', @ischar)
addOptional(inPars,'gpu',0, @isnumeric)
addOptional(inPars,'sort','r2s', @ischar)
parse(inPars,images,varargin{:})

rhdr=images.hdr;
cdata=images.cdata;
EchoSpacing=images.EchoSpacing;

if inPars.Results.ROI==0
    disp('Select ROI for processing')
    imagesc(abs(cdata(:,:,1,1,1)))
    roi=round(getrect);
    close 
else
    roi=inPars.Results.ROI;
end

if inPars.Results.driftROI==1
    disp('Select ROI to calculate drift')
    imagesc(abs(cdata(:,:,1,1,1)))
    droi=round(getrect);
    close 
else
    droi=inPars.Results.driftROI;
end

modelflag=inPars.Results.model;
gpuflag  =inPars.Results.gpu;
firstImage=inPars.Results.firstImage;
numPeaks=inPars.Results.numPeaks;
sortflag=inPars.Results.sort;

alpha=-inPars.Results.alpha; % *fudged sign correction- should be corrected elsewhere in the script

threshold=inPars.Results.thresh; %arbitrary threshold for processing

workers=inPars.Results.procs;
poolopenflag=matlabpool('size'); %check to see if matlabpool is open
if workers>0 && poolopenflag==0 % If you want to use multiple procs and the pool isnt open, open it.
    matlabpool open
end

%% Read header for expected array dimensions, preallocate space for arrays
%dicomdict('set','/FUS/matlab/gems-dicom-dictHDX.txt');
%NumTimePoints=rhdr(1,1,1).NumberOfTemporalPositions;
NumTimePoints=numel(images.cdata(1,1,1,1,:));

%NumTimePoints=rhdr(1,1,1).NumberOfAcquisitions;
MatSizeX=numel(cdata(1,:,1,1,1));
MatSizeY=numel(cdata(:,1,1,1,1));
NumEchoes=numel(images.EchoTimes);%rhdr(1,1,1).EchoTrainLength;
%NumSlices=rhdr(1,1,1).Private_0021_104f;
NumSlices=rhdr(1,1,1).LocationsInAcquisition;

f0=rhdr.ImagingFrequency;
% execute pixel by pixel CSI processing (uber big for loops...) 
t2star=zeros(MatSizeY,MatSizeX,NumSlices,NumTimePoints,numPeaks);
ppm=t2star;
amp=ppm;
t1=amp;
tcdata=ones([size(cdata),numPeaks]);
deltaT=zeros(MatSizeY,MatSizeX,NumSlices,NumTimePoints);

%% Reshape for parfor, perform ARMA processing
h1=waitbar(0,'CSI Processing');
cdata2=cdata(roi(2):roi(2)+roi(4),roi(1):roi(1)+roi(3),:,:,:);
nrows=numel(cdata2(:,1,1,1,1));
ncols=numel(cdata2(1,:,1,1,1));
cdata2=reshape(cdata2,(nrows*ncols),NumSlices,NumEchoes,NumTimePoints);

ppmtmp2=zeros(nrows*ncols,NumSlices,NumTimePoints,numPeaks);
t2startmp2=ppmtmp2;
amptmp2=ppmtmp2;
t1tmp2=ppmtmp2;

switch 1
    case gpuflag==1
        
        cdata2=reshape(permute(cdata2,[3 1 2]),1,[]);
        
        
        % inialize device data
        %d_fid     = gpuArray( complex(cdata2,cdata2) );
        d_fid     = gpuArray( cdata2     );
        d_ppm     = gpuArray( ppmtmp2    );
        d_t2star  = gpuArray( t2startmp2 );
        d_mag     = gpuArray( amptmp2    );
        d_phase   = gpuArray( amptmp2    );
        
        PronyKernel = parallel.gpu.CUDAKernel('PronyKernel.ptx', 'PronyKernel.cu');
        
        threadsPerBlock = 256;
        PronyKernel.ThreadBlockSize=[threadsPerBlock  1];
        npixel = (nrows*ncols);
        blocksPerGrid = (npixel  + threadsPerBlock - 1) / threadsPerBlock;
        PronyKernel.GridSize=[ceil(blocksPerGrid)  1];
        
        % echo info
        disp( sprintf('NumTimePoints  %d NumSlices %d NumEchoes %d numPeaks %d npixel %d',NumTimePoints,NumSlices,NumEchoes,numPeaks,npixel) );
        
        [d_ppm, d_t2star, d_mag, d_phase ] = feval(PronyKernel,real(d_fid),imag(d_fid), d_ppm, d_t2star, d_mag,d_phase   ,EchoSpacing,f0,NumEchoes,numPeaks,npixel);
        
        ppmtmp2=permute(reshape(1000*gather(d_ppm),[numPeaks,ncols,nrows]),[2 3 1]);
        t2startmp2=permute(reshape(gather(d_t2star)/1000,[numPeaks,ncols,nrows]),[2 3 1]);
        
        magtmp=gather(d_mag);
        phtmp=gather(d_phase);
        amptmp2=permute(reshape(complex(magtmp.*cos(phtmp),magtmp.*sin(phtmp)),[numPeaks,ncols,nrows]),[2 3 1]);
        
    case workers>0
        armafun=@(xx) ooCsi_arma2(xx,EchoSpacing,f0,numPeaks,modelflag);
        for tt=1:NumTimePoints;
            for ss=1:NumSlices;
                parfor (ll=1:nrows*ncols,workers)
                    %for ll=1:nrows*ncols
                    if abs(cdata2(ll,ss,1,tt))>threshold
                        fid=squeeze(cdata2(ll,ss,:,tt));
                        warning('off','MATLAB:toeplitz:DiagonalConflict');
                        warning('off','MATLAB:legend:IgnoringExtraEntries');
                        warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
                        [ppmtmp,t2stmp,amptmp,~]=armafun(fid);
                        ppmtmp2(ll,ss,tt,:)=ppmtmp;
                        t2startmp2(ll,ss,tt,:)=t2stmp;
                        amptmp2(ll,ss,tt,:)=amptmp;
                    else
                        t2startmp2(ll,ss,tt,:)=ones(numPeaks,1)*NaN;
                        ppmtmp2(ll,ss,tt,:)=ones(numPeaks,1)*NaN;
                        amptmp2(ll,ss,tt,:)=ones(numPeaks,1)*NaN;
                    end
                end
                
            end
            waitbar((double(tt*ss))/(double(NumSlices*NumTimePoints)),h1,'CSI Processing');
        end
        
        delete(h1)
        
    case workers==0
        armafun=@(xx) ooCsi_arma(xx,EchoSpacing,f0,numPeaks,modelflag);
        for tt=1:NumTimePoints;
            for ss=1:NumSlices;
                for ll=1:nrows*ncols
                    if abs(cdata2(ll,ss,1,tt))>threshold
                        fid=squeeze(cdata2(ll,ss,:,tt));
                        warning('off','MATLAB:toeplitz:DiagonalConflict');
                        warning('off','MATLAB:legend:IgnoringExtraEntries');
                        warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
                        [ppmtmp,t2stmp,amptmp,~]=armafun(fid);
                        ppmtmp2(ll,ss,tt,:)=ppmtmp;
                        t2startmp2(ll,ss,tt,:)=t2stmp;
                        amptmp2(ll,ss,tt,:)=amptmp;
                    else
                        t2startmp2(ll,ss,tt,:)=ones(numPeaks,1)*NaN;
                        ppmtmp2(ll,ss,tt,:)=ones(numPeaks,1)*NaN;
                        amptmp2(ll,ss,tt,:)=ones(numPeaks,1)*NaN;
                    end
                end
                
            end
            waitbar((double(tt*ss))/(double(NumSlices*NumTimePoints)),h1,'CSI Processing');
        end
        delete(h1)
end

%%
ppmtmp3=reshape(ppmtmp2,[nrows,ncols,NumSlices,NumTimePoints,numPeaks]);
t2startmp3=reshape(t2startmp2,[nrows,ncols,NumSlices,NumTimePoints,numPeaks]);
amptmp3=reshape(amptmp2,[nrows,ncols,NumSlices,NumTimePoints,numPeaks]);
t1tmp3=reshape(t1tmp2,[nrows,ncols,NumSlices,NumTimePoints,numPeaks]);
sigtmp3=zeros([size(ppm),NumEchoes]);

ppm(roi(2):roi(2)+roi(4),roi(1):roi(1)+roi(3),:,:,:)=ppmtmp3;
t2star(roi(2):roi(2)+roi(4),roi(1):roi(1)+roi(3),:,:,:)=t2startmp3;
amp(roi(2):roi(2)+roi(4),roi(1):roi(1)+roi(3),:,:,:)=amptmp3;
t1(roi(2):roi(2)+roi(4),roi(1):roi(1)+roi(3),:,:,:)=t1tmp3;
%%
msk=permute(abs(cdata(:,:,:,1,:))>threshold,[1 2 3 5 4]);
msk=repmat(msk,[1 1 1 1 numPeaks]);

%calculate tmap
%Decide which image to use as baseline %Use first peak by default
if NumTimePoints>=(firstImage+1);
    for mm=1:NumSlices
        for nn=firstImage:NumTimePoints
            deltat(:,:,mm,nn,1)=(ppm(:,:,mm,nn,1)-ppm(:,:,mm,firstImage,1))/alpha;
        end
        if droi ~= 0
            drift(1,1,mm,:,1)=squeeze(mean(mean(imcrop2(deltat(:,:,mm,:,1),droi),2),1));
            drift2=repmat(drift,numel(deltat(:,1,1,1,1)),numel(deltat(1,:,1,1,1)));
        else
            drift2=zeros(size(deltat));
        end
    end
    data.deltaT=(deltat-drift2).*msk(:,:,:,:,1);
end

if droi~=0
    data.driftROI=droi;
    data.drift=drift;
end


%% Sort data

switch 1
    case strcmp(sortflag,'amp')||strcmp(sortflag,'+amp')
        sortarray=nan2zero(amp);
        sortstring='ascend';
    case strcmp(sortflag,'-amp')
        sortarray=nan2zero(amp);
        sortstring='descend';
    case strcmp(sortflag,'ppm')||strcmp(sortflag,'+ppm')
        sortarray=nan2zero(ppm);
        sortstring='ascend';
    case strcmp(sortflag,'-ppm')
        sortstring='descend';
        sortarray=nan2zero(ppm);
    case strcmp(sortflag,'r2s')||strcmp(sortflag,'+r2s')
        sortstring='descend';
        sortarray=nan2zero(t2star);
    case strcmp(sortflag,'-r2s')
        sortarray=nan2zero(t2star);
        sortstring='ascend';
    case strcmp(sortflag,'t2s')||strcmp(sortflag,'+t2s')
        sortarray=nan2zero(t2star);
        sortstring='ascend';
    case strcmp(sortflag,'-t2s')
        sortarray=nan2zero(t2star);
        sortstring='descend';
end

[~,sortindx]=sortlind(sortarray,5,sortstring);

data.cdata=cdata;
data.t2star=1000*t2star(sortindx).*msk;
data.ppm=ppm(sortindx).*msk;
data.amp=amp(sortindx).*msk;
data.t1=t1(sortindx).*msk;
% Need to output signal
%data.tcdata=tcdata;
data.ROI=roi;
end

