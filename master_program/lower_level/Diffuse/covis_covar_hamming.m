function [cov,E1,E2,rc] = covis_covar_hamming(a,b,r,window,overlap)

% 29 Aug. 2018, drj12@uw.edu, Error in earlier version corrected
% Adapted from covis_cor_win.m. This version does not normalize by energy
% and applies a Hamming window. This function returns separate energies for
% the two signals.

%   Fixed window size covariance of two complex series
%
%   INPUTS:
%      a,b = complex series
%      r = range vector
%      window = window size (number of samples)
%      overlap = window overlap size (number of samples)
%
%   outputs:
%      rc = bin range
%      cov = covariance
%      E1 = energy of windowed a
%      E2 = Energy of windowed b
%
%------


[N,M]=size(a);
ham_win = hamming(window)*ones(1,M)/sum(hamming(window));
nbins=floor(N/(window-overlap))-1; % Last range bin is removed
cab=a.*conj(b);
cov=zeros(nbins,M);
E1=zeros(nbins,M);
E2=zeros(nbins,M);
rc=zeros(nbins,1);
for m=1:nbins
    j=(m-1)*(window-overlap)+1;
    k=j+window-1;
    if (k>N) k=N; end
    i=j:k;
    cabjk = cab(i,:);
    cov(m,:)=sum(cabjk.*ham_win);
    E1(m,:) = sum(abs(a(i,:)).^2.*ham_win);
    E2(m,:) = sum(abs(b(i,:)).^2.*ham_win);
    rc(m)=mean(r(i));
end

end




