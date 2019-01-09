function [num,den]=OOcsi_arma(fid,esp,nPeaks,model)
nPeaks0=nPeaks;
esp=10^-3*esp;
if model=='stmcb'
    coeffun=@stmcb;
elseif model=='prony'
    coeffun=@prony;
end
r2s(1)=-1;
mod_index=0;

while min(r2s)<0
    
    % Reduce model order if negative t2star/r2star
    npeaks=nPeaks0;% npeaks=nPeaks0-mod_index;
    
    % Find coefficients/roots
    [num,den]=coeffun(fid,npeaks-1,npeaks);
    mod_index=mod_index+1;
end