function [covar,rc]=autocovar_win(z,r,window,overlap,lag)
%
%   Average over a window of the unit-lag
%   autocovariance of a complex series.
%
%   INPUTS:
%      z = complex series
%      r = range vector
%      window = window (number of samples)
%      overlap = window overlap (number of samples)
%      lag = ping-lag used to calculate covariance
%
%   OUTPUT:
%	covar = covariance
%	rc = positions of covariance values (m)

[N,M]=size(z);
nbins=floor(N/(window-overlap));
covar=zeros(nbins,M);
rc=zeros(nbins,1);
for m=1:nbins
    j=(m-1)*(window-overlap)+1;
    k=j+window-1;
    
    if j < 1
        j = 1;
    end
    if k>N-lag
        k=N-lag;
    end
    if j > k
        j = k;
    end
    
    i=j:k;
    covar(m,:)=mean(z(i,:).*conj(z(i+lag,:)));
    rc(m)=mean(r(i));
    
end

end