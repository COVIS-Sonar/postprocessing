% This function calculates the volume of the plume above Inferno by
% counting the pixels with VSS within a seleted range

function [vol,vss] = get_plume_volume_vss_fun(grid,center_coor,zmin,zmax,vss_low,vss_high)
% Input:
% covis: Matlab structure that includes gridded Imaging data
% xs: x-coordinate of sourve vent (e.g., Inferno) (m)
% ys: y-coordinate of source vent (m)
% Ro: search radius for plume signals (m)
% vss_low: low end of the vss range used to define plume signals (dB)
% vss_high: high end of the vss range used to deefine plume signals (dB)

% Output:
% vol: plume volume (m^3)
% vss: volume averaged vss (dB)


% Example
% xs = -13.27;
% ys = -0.77;
% Ro = 4;
% vss_low = -60;
% vss_high = -20;

% xrange = [xs-Ro,xs+Ro];
% yrange = [ys-Ro,ys+Ro];
% zrange = [zmin, zmax];

alpha = 0.1;
z_vent = 4;
b0 = 0.5;

Ig = grid.Id;
Ig(Ig==10^-9) = nan;
vss_in = 10*log10(Ig);
xc = center_coor(:,1);
yc = center_coor(:,2);
zc = center_coor(:,3);

xg=grid.x;
yg=grid.y;
zg=grid.z;
xg2 = squeeze(xg(:,:,1));
yg2 = squeeze(yg(:,:,1));
zg1 = squeeze(zg(1,1,:));
zind = zg1>=zmin&zg1<=zmax;
zg1_sub = zg1(zind);
vss_sub = vss_in(:,:,zind);

b = (zg1_sub - z_vent)*6/5*alpha*2+b0;

dz=grid.spacing.dz;
dx=grid.spacing.dx;
dy=grid.spacing.dy;

vol_temp = zeros(1,length(zg1_sub));
vss_temp = zeros(1,length(zg1_sub));
for k = 1:length(zg1_sub)
    vss1 = squeeze(vss_sub(:,:,k));
    xc1 = interp1(zc,xc,zg1_sub(k));
    yc1 = interp1(zc,yc,zg1_sub(k));
    b1 = b(k);
    xrange = [xc1-b1,xc1+b1];
    yrange = [yc1-b1,yc1+b1];
    cond1 = xg2>=xrange(1)&xg2<=xrange(2);
    cond2 = yg2>=yrange(1)&yg2<=yrange(2);
    cond3 = vss1>=vss_low&vss1<=vss_high;
    ii = find(cond1&cond2&cond3);
    vol_temp(k) = length(ii)*dx*dy*dz;
    vss_temp(k) = nanmean(vss1(ii));
end
vol = nansum(vol_temp);
vss = nanmean(vss_temp);

% cond1 = xg>=xrange(1)&xg<=xrange(2);
% cond2 = yg>=yrange(1)&yg<=yrange(2);
% cond3 = zg>=zrange(1)&zg<=zrange(2);
% cond4 = vss_in>=vss_low&vss_in<=vss_high;

% ii = find(cond1&cond2&cond3&cond4);
% vol = length(ii)*dx*dy*dz;
% vss = nanmean(vss_in(ii));





