function grd_out = l3grid_doppler(grd_in,grd_out)
% Grids the values in 2-D grid 'grd_in' with coordinates (x,y,z) onto a uniform 
% 3-D grid 'grd_out' with coordinates (xg,yg,zg) using nearest-neighbor linear interpolation. 
% Data is interpolated onto the neighboring grid cell only. 
% The positions (x,y,z) may be non-uniform and non-monotinic data.  
% Interpolation is done using a linear weighting matrix.
% An updated weighting matrix is returned after each call, 
% so the function can be called multiple time on the same grid, 
% adding new data to the grid on each call. 
% After the final call, the gidded data should be normalized 
% by the weighting matrix, like this
%   n=find(wg); Ig(n)=Ig(n)./wg(n);
%
% ---
% Version 1.0 - cjones@apl.washington.edu 06/2010
% Version 1.1 - xupeng_66@hotmail.com 05/2011
% Version 1.2 - xupeng_66@hotmail.com 06/2011 (adding a grid for the radial
% velocity vr)
% Version 1.3 - guangyux@uw.edu 09/2018 (changed input and output structures)

x = grd_in.x;
y = grd_in.y;
z = grd_in.z;
xg = grd_out.x;
yg = grd_out.y;
zg = grd_out.z;

[Ny,Nx,Nz]=size(xg);

dy = abs(yg(2,1,1)-yg(1,1,1));
dx = abs(xg(1,2,1)-xg(1,1,1));
dz = abs(zg(1,1,2)-zg(1,1,1));

xg=xg(:);
yg=yg(:);
zg=zg(:);

xmin=min(xg);
ymin=min(yg);
zmin=min(zg);
xmax=max(xg);
ymax=max(yg);
zmax=max(zg);

% use only finite points within the grid
if strcmp(grd_out.type,'doppler velocity');
    vr_cov = grd_in.vr_cov;
    vr_vel = grd_in.vr_vel;
    vz_cov = grd_in.vz_cov;
    std = grd_in.std;
    covar = grd_in.covar;
    ii=find(isfinite(vr_cov)&(x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)&(z>=zmin)&(z<=zmax));
    if (~isempty(ii))
        x=x(ii);
        y=y(ii);
        z=z(ii);
        vr_cov=vr_cov(ii);
        vr_vel=vr_vel(ii);
        vz_cov=vz_cov(ii);
        std = std(ii);
        covar = covar(ii);
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
                    grd_out.vr_cov(p)=grd_out.vr_cov(p)+w.*vr_cov(:); % covariance-averaged line-of-sight velocity
                    grd_out.vr_vel(p)=grd_out.vr_vel(p)+w.*vr_vel(:); % velocity-averaged line-of-sight velocity 
                    grd_out.vz_cov(p)=grd_out.vz_cov(p)+w.*vz_cov(:); % pseudo vertical velocity 
                    grd_out.std(p)=grd_out.std(p)+w.*std(:); % velocity standard deviation
                    grd_out.covar(p) = grd_out.covar(p)+w.*covar(:); % covariance function
                    grd_out.w(p)=grd_out.w(p)+w; % weight function
                end
            end
        end
    end
elseif strcmp(grd_out.type,'intensity')
    I = grd_in.I;
    std = grd_in.std;
    ii=find(isfinite(I)&(x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)&(z>=zmin)&(z<=zmax));
    if (~isempty(ii))
        x=x(ii);
        y=y(ii);
        z=z(ii);
        I=I(ii);
        std = std(ii);
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
                    grd_out.I(p)=grd_out.I(p)+w.*I(:); % ping-averaged volume backscattering coefficient
                    grd_out.std(p)=grd_out.std(p)+w.*std(:); % standard deviation of volume backscattering coefficient
                    grd_out.w(p)=grd_out.w(p)+w; % weight function
                end
            end
        end
    end
else
    error('No Doppler grid is found')
end
end






