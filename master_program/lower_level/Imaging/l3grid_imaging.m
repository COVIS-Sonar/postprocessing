function grd_out = l3grid_imaging(grd_in,grd_out)


% Grids the values in grid 'grd_in' with coordinates (x,y,z) onto a uniform 
% 3-D grid grd_out with coordinates(xg,yg,zg) using weighted nearest-neighbor 
% linear interpolation. Data is interpolated onto the neighboring grid cell
% only. The positions (x,y,z) may be non-uniform and non-monotinic data.  
% Interpolation is done using a linear weighting matrix.
% An updated weighting matrix is returned after each call, 
% so the function can be called multiple time on the same grid, 
% adding new data to the grid on each call. 
% After the final call, the gidded data should be normalized 
% by the weighting matrix, like this
%   n=find(wg); Ig(n)=Ig(n)./wg(n);
%
% ---
% Version 1.0 by guangyux@uw.edu (Nov 13,2019)
%  --based on the original code written by Chris Jones in 2010

x = grd_in.x;
y = grd_in.y;
z = grd_in.z;
xg1 = grd_out.x;
yg1 = grd_out.y;
zg1 = grd_out.z;

[Ny,Nx,Nz]=size(xg1);


dx = abs(xg1(1,2,1)-xg1(1,1,1));
dy = abs(yg1(2,1,1)-yg1(1,1,1));
dz = abs(zg1(1,1,2)-zg1(1,1,1));

% The following three lines are added to ensure correct indexing in the
% calculation of interpolation weight
xg=squeeze(xg1(1,:,1));
yg=squeeze(yg1(:,1,1));
zg=squeeze(zg1(1,1,:));

xg = xg(:);
yg = yg(:);
zg = zg(:);

xmin=min(xg);
ymin=min(yg);
zmin=min(zg);
xmax=max(xg);
ymax=max(yg);
zmax=max(zg);

Ia = grd_in.Ia;
Id = grd_in.Id;
Ia_filt = grd_in.Ia_filt;
Id_filt = grd_in.Id_filt;
Kp = grd_in.Kp;
% use only finite points within the grid
ii=find(isfinite(Id)&isfinite(Ia)&(x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)&(z>=zmin)&(z<=zmax));

if(~isempty(ii))

x=x(ii);
y=y(ii);
z=z(ii);
Ia=Ia(ii);
Id=Id(ii);
Ia_filt = Ia_filt(ii);
Id_filt = Id_filt(ii);
Kp = Kp(ii);

i(:,1)=floor((x(:)-xmin)/dx)+1;
i(:,2)=ceil((x(:)-xmin)/dx)+1;
j(:,1)=floor((y(:)-ymin)/dy)+1;
j(:,2)=ceil((y(:)-ymin)/dy)+1;
k(:,1)=floor((z(:)-zmin)/dz)+1;
k(:,2)=ceil((z(:)-zmin)/dz)+1;

for n=1:2
   wx=1-abs(xg(i(:,n))-x(:))/dx;
   for m=1:2
      wy=1-abs(yg(j(:,m))-y(:))/dy;
      for l=1:2
         wz=1-abs(zg(k(:,l))-z(:))/dz;
         w=sqrt(wx.^2+wy.^2+wz.^2);
         p=sub2ind([Ny,Nx,Nz],j(:,m),i(:,n),k(:,l));
         grd_out.Ia(p)=grd_out.Ia(p)+w.*Ia(:); % averaged backscatter intensity
         grd_out.Id(p)=grd_out.Id(p)+w.*Id(:); % differenced backscatter intensity
         grd_out.Ia_filt(p)=grd_out.Ia_filt(p)+w.*Ia_filt(:); % filtered averaged backscatter intensity
         grd_out.Id_filt(p)=grd_out.Id_filt(p)+w.*Id_filt(:); % filtered differenced backscatter intensity
         grd_out.Kp(p) = grd_out.Kp(p)+w.*Kp(:);
         grd_out.w(p)=grd_out.w(p)+w; % weight function
      end
   end
end

end
