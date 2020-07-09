% This function is used to create 3D acoustic images of volume scattering 
% strength (VSS) from the gridded Imaging-mode data

% version 1.0 by guangyux@uw.edu (Oct 19, 2019)
%  --based on the original code written by Chris Jones in 2010
% version 2.0 by guangyux@uw.edu (Feb 20, 2020)
%  --new color schemes for bathymetry and VSS isosurfaces
%  --new VSS isosurface values
%  --new angles of viewpoint

function covis_imaging_plot(covis)
% Input:
%  covis: the Matlab structure that contains the gridded Imaging data and
%         metadata
% Output:
%  a 3D acoustic plume image


% set up constant parameters
threshold = 3; % threshold for differentiating plumes and structures (dB)


% load data and grid
swp_name = covis.sweep.name;
swp_date = datenum(swp_name(7:21),'yyyymmddTHHMMSS');
xg = covis.grid.x;
yg = covis.grid.y;
zg = covis.grid.z;
vg_a = covis.grid.Ia;
vg_d = covis.grid.Id;
dvg = vg_d./vg_a;
vg = vg_d;
vg(10*log10(dvg)<threshold) = 10^-9;




% load bathymetry data
if swp_date<=datenum(2019,7,6)
    bathy_name = sprintf('covis_bathy_2018.mat');
elseif swp_date<=datenum(2019,11,23)
    bathy_name = sprintf('covis_bathy_2019a.mat');
else
    bathy_name = sprintf('covis_bathy_2019b.mat');
end
bathy_file = find_input_file(bathy_name);
covis_bathy = load(bathy_file);
xb = covis_bathy.covis.grid.x;
yb = covis_bathy.covis.grid.y;
zb = covis_bathy.covis.grid.v;
rb = sqrt(xb.^2+yb.^2);
zb(rb<4) = nan;



% plot 3D image
figure
% add bathy
pbathy=surf(xb,yb,zb);
axis image;
caxis([-2 4]);
caxis('manual');
pbathy.FaceColor='interp';
pbathy.EdgeColor='none';
pbathy.DiffuseStrength=0.9;
pbathy.BackFaceLighting='lit';
pbathy.SpecularColorReflectance=0;
pbathy.SpecularExponent=20;
%pbathy.SpecularStrength=0.3;
pbathy.SpecularStrength=0.9;
temp=pink(16);
mymap=flipud(temp(6:13,:));
colormap(mymap)
hold on
contour3(xb,yb,zb,[-2:0.2:4],'LineColor',[0.4 0.4 0.4])   
isosurf{1}.color = [77,153,204]/255;
isosurf{2}.color = [128,77,128]/255;
isosurf{3}.color = [153,5,13]/255;
isosurf{1}.value = -60;
isosurf{2}.value = -50;
isosurf{3}.value = -40;
isosurf{1}.alpha = 0.1;
isosurf{2}.alpha = 0.2;
isosurf{3}.alpha = 0.3;
for n=1:length(isosurf)
    v = vg;
    eps = nan;
    m = find(v==0);
    v(m) = eps;  % remove zeros
    v = 10*log10(v);
    surf_value = isosurf{n}.value;
    surf_color = isosurf{n}.color;
    surf_alpha = isosurf{n}.alpha;
    p = patch(isosurface(xg, yg, zg, v, surf_value));
    isonormals(xg, yg, zg, v, p);
    set(p,'EdgeColor','none','FaceColor',surf_color,'FaceAlpha',surf_alpha);
end
daspect([1 1 1])
xlim([-20 1]);
ylim([-8 12]);
zlim([-2 15]);
plot3([0,0],[0,0],[0,4.2],'y','linewidth',4);
hold off;
az = 76.5;
el = 26;
view([az,el])
camlight('headlight');
xlabel('Easting of COVIS ( m )')
ylabel('Northing of COVIS ( m )')
zlabel('Height above COVIS base ( m )')
title(covis.sweep.name);
h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
end

