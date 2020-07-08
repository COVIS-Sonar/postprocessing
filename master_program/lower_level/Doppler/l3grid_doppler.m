function grd_out = l3grid_doppler(grd_in,grd_out)
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
if strcmp(grd_out.type,'doppler velocity')
    vr_cov = grd_in.vr_cov;
    vr_vel = grd_in.vr_vel;
    std = grd_in.std;
    covar = grd_in.covar;
    vr_cov(vr_cov==0) = nan;
    vr_vel(vr_vel==0) = nan;
    std(std==0) = nan;
    covar(covar==0) = nan;
    ii=find(isfinite(vr_cov)&(x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)&(z>=zmin)&(z<=zmax));
    if (~isempty(ii))
        x=x(ii);
        y=y(ii);
        z=z(ii);
        vr_cov=vr_cov(ii);
        vr_vel=vr_vel(ii);
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
                    grd_out.std(p)=grd_out.std(p)+w.*std(:); % velocity standard deviation
                    grd_out.covar(p) = grd_out.covar(p)+w.*covar(:); % covariance function
                    grd_out.w(p)=grd_out.w(p)+w; % weight function
                end
            end
        end
    end
    grd_out.vr_cov(isnan(grd_out.vr_cov)) = 0;
    grd_out.vr_vel(isnan(grd_out.vr_vel)) = 0;
    grd_out.std(isnan(grd_out.std)) = 0;
    grd_out.covar(isnan(grd_out.covar)) = 0;
elseif strcmp(grd_out.type,'intensity') || strcmp(grd_out.type,'intensity_win')
    Id = grd_in.Id;
    Id_filt = grd_in.Id_filt;
    Id(Id==10^-9) = nan;
    Id_filt(Id==10^-9) = nan;
    ii=find(isfinite(Id)&(x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)&(z>=zmin)&(z<=zmax));
    if (~isempty(ii))
        x=x(ii);
        y=y(ii);
        z=z(ii);
        Id=Id(ii);
        Id_filt = Id_filt(ii);
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
                    grd_out.Id(p)=grd_out.Id(p)+w.*Id(:); % ping-averaged volume backscattering coefficient
                    grd_out.Id_filt(p)=grd_out.Id_filt(p)+w.*Id_filt(:); % standard deviation of volume backscattering coefficient
                    grd_out.w(p)=grd_out.w(p)+w; % weight function
                end
            end
        end
    end
    grd_out.Id(isnan(grd_out.Id)) = 10^-9;
    grd_out.Id_filt(isnan(grd_out.Id_filt)) = 10^-9;
else
    error('No Doppler grid is found')
end
end






