function [num,den,rtz]=csi_armaMGE(fid,order,model)
nPeaks0=order;
if model=='stmcb'
    coeffun=@stmcb;
elseif model=='prony'
    coeffun=@prony;
end
r2s(1)=-1;
mod_index=0;

% num=nan*ones(1,nPeaks0);
% den=nan*ones(1,nPeaks0+1);
% rtz=nan*ones(1,nPeaks0);

while min(r2s)<0 && (nPeaks0-mod_index)>=1
    
    num=zeros(1,nPeaks0);
    den=zeros(1,nPeaks0+1);
    rtz=zeros(1,nPeaks0);
    
    % Reduce model order if negative t2star/r2star
    npeaks=nPeaks0-mod_index;
    
    % Find coefficients/roots
    warning('off','MATLAB:toeplitz:DiagonalConflict');
    warning('off','MATLAB:legend:IgnoringExtraEntries');
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
    [num(1:npeaks),den(1:npeaks+1)]=coeffun(fid,[1 zeros(1,numel(fid)-1)],npeaks-1,npeaks);
    
    %num=flipdim(num,2);
    %den=flipdim(den,2);
    %r2s=real(log(rtz));
    %rtz(1:npeaks)=transpose(roots(den(1:npeaks+1)));
    
    rtz(1:npeaks)=roots(den(1:npeaks+1));
    r2s=-real(log(rtz));
    mod_index=mod_index+1;
end

end