function y=gaussian2D(beta,x)
% returns value of gaussian for beta and x as defined by
%		y=beta(1)*exp(-x^2/beta(2)^2)+beta(3)
beta(2)=1e6*beta(2);
beta(4)=1e6*beta(4);
beta(5)=1e6*beta(5);
r2n=(x(:,1)-beta(4)).^2 + (x(:,2)-beta(5)).^2;
y=beta(1).*exp(-r2n./(beta(2).^2))+beta(3);