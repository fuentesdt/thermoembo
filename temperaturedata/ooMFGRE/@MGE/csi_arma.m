function [num,den]=csi_armaMGE(fid,order,model)
nPeaks0=obj.nPeaks;
if obj.model=='stmcb'
    coeffun=@stmcb;
elseif obj.model=='prony'
    coeffun=@prony;
end
%r2s(1)=-1;
%mod_index=0;

%while min(r2s)<0
    
    % Reduce model order if negative t2star/r2star
    npeaks=nPeaks0;% npeaks=nPeaks0-mod_index;
    
    % Find coefficients/roots
    [num,den]=coeffun(fid,npeaks-1,npeaks);
    %mod_index=mod_index+1;
%end
end