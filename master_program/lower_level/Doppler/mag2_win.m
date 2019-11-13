function [zm,rc]=mag2_win(z,r,window,overlap)
%
%   Average over a window of the magnitude squared of a complex series
%
%   INPUTS:
%      z = complex series
%      r = range vector
%      window = window (number of samples)
%      overlap = window overlap (number of samples)
%   option = 'unbiased', 'biased', or 'coeff'
%
%   OUTPUT:
%	zm = window averaged magnitude squared of z
%	rc = positions of wabs values (m)

[N,M]=size(z);
nbins=floor(N/(window-overlap));
zm=zeros(nbins,M);
rc=zeros(nbins,1);
for m=1:nbins
    j=(m-1)*(window-overlap)+1;
    k=j+window-1;
    
    if j < 1
        j = 1;
    end
    if k>N
        k=N;
    end
    if j > k
        j = k;
    end
    
    i=j:k;
    zm(m,:)=mean(abs(z(i,:).^2));
    if k <= N
        rc(m)=mean(r(i));
    else
        rc(m)=r(N);
    end
    
end
end