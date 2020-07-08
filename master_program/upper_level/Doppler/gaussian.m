function y=gaussian(beta,x,fac)
% returns value of gaussian for beta and x as defined by
%		y=beta(1)*exp(-(x-beta(4))^2/beta(2)^2)+beta(3)
%   which is a form of y=A*exp(-r^2/b^2)+c
%           beta(1) = A
%           beta(2) = b
%           beta(3) = c  (background)
%           beta(4) = x0 (centerline location)
%
beta(2)=beta(2)/fac;
beta(4)=beta(4)/fac;
y=beta(1).*exp(-((x-beta(4)).^2)./(beta(2).^2))+beta(3);