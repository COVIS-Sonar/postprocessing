function y=gaussian(beta,x)
% returns value of gaussian for beta and x as defined by
%		y=beta(1)*exp(-(x-beta(4))^2/beta(2)^2)+beta(3)
%   which is a form of y=A*exp(-r^2/b^2)+c
%           beta(1) = A
%           beta(2) = b
%           beta(3) = c  (background)
%           beta(4) = x0 (centerline location)
%
beta(2)=1e6*beta(2);
beta(4)=1e6*beta(4);
y=beta(1).*exp(-((x-beta(4)).^2)./(beta(2).^2))+beta(3);