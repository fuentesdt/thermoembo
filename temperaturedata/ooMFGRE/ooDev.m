clear all
profile on
tic
tmp=MGE('/mnt/FUS4/data2/CJM/MSC_SPIO/20140304/Raw/tess7T.oC2/4/');
tmp.processingROI=[24 29 110 42];
%tmp.setProcessingROI
tmp.numProcs=12;
tmp.modelOrder=1;
tmp.calcArmaCoeffs;
imagesc(abs(tmp.amp(:,:,12,1,1)));
imagesc(tmp.t2star(:,:,12,1,1));
imagesc(tmp.ppm(:,:,12,1,1));
time1=toc;
profile off
profile viewer
whos
%%
clear all
tic
tmp2=importMFGRE('/mnt/FUS4/data2/CJM/MSC_SPIO/20140304/Raw/tess7T.oC2/4/');
tmp3=arma(tmp2,'ROI',[24 29 110 42],'procs',12);
imagesc(abs(tmp3.amp(:,:,12,1,1)));
imagesc(tmp3.t2star(:,:,12,1,1));
imagesc(tmp3.ppm(:,:,12,1,1));
time2=toc;
whos
% a=tmp.armaNumerator;
% b=tmp.armaDenominator;
% 
% 
% tic
% tmp.calcArmaCoeffs;
% toc
% 
% a2=tmp.armaNumerator;
% b2=tmp.armaDenominator;
% 
% isequal(a,a2)
% isequal(b,b2)
% %mfgreDCM='/FUS4/data2/CJM/MnCuS/20131024/e162/s1136';
% %mfgreBruk3D='/mnt/FUS4/data2/CJM/MSC_SPIO/20131213/Raw/tess7T.nf1/5/';

