% This function reads gridded Imaging data to plot 3D plume images

% version 1.0 by guangyux@uw.edu (Oct 19, 2019)
%  --based on the original code written by Chris Jones in 2010


function covis_imaging_plot(covis)
% Input:
% covis: structure array containing the gridded Imaging data and associated
% meta data

swp_name = covis.sweep.name;
year = swp_name(7:10);
xg = covis.grid.x;
yg = covis.grid.y;
zg = covis.grid.z;
vg = covis.grid.Id_filt;

bathy_name = sprintf('covis_bathy_%s.mat',year);
bathy = load(bathy_name);
covis_bathy = bathy.covis;
xb = covis_bathy.grid.x;
yb = covis_bathy.grid.y;
zb = covis_bathy.grid.v;

% plot 3D image
figure
surf(xb,yb,zb)
shading interp
colormap('summer');
axis image;
hold on;
caxis([-1 3]);
caxis('manual');
contour3(xb,yb,zb,20,'k');
isosurf{1}.color = [0,0,1];
isosurf{2}.color = [0.5 0 0.5];
isosurf{3}.color = [1,0,0];
isosurf{1}.value = -60;
isosurf{2}.value = -50;
isosurf{3}.value = -40;
isosurf{1}.alpha = 0.2;
isosurf{2}.alpha = 0.4;
isosurf{3}.alpha = 0.6;
%for n = 3
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
zlim([-2 15]);
plot3([0,0],[0,0],[0,4.2],'y','linewidth',4);
hold off;
az = 76.5;
el = 26;
view([az,el])
xlabel('Easting of COVIS ( m )')
ylabel('Northing of COVIS ( m )')
zlabel('Height above COVIS feet ( m )')
end
