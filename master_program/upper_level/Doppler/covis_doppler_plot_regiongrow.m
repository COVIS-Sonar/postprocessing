% This function is used to create 3D acoustic images of volume scattering 
% strength (VSS) from the gridded Imaging-mode data


function covis_doppler_plot_regiongrow(covis)
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
xg = covis.grid{3}.x;
yg = covis.grid{3}.y;
zg = covis.grid{3}.z;
xgg = squeeze(xg(:,:,1));
ygg = squeeze(yg(:,:,1));
zgg = squeeze(zg(1,1,:));
vg_a = covis.grid{3}.Ia;
vg_d = covis.grid{3}.Id;
dvg = vg_d./vg_a;
vg = vg_d;
vg(10*log10(dvg)<threshold) = 10^-9;


% smooth data 
vg_plot = smooth3(vg,'gaussian',7,0.3);


% use the method region growing to mask out noise and artifacts
R = 4;
vg_thresh = exp(-8);
x0 = -13.27;
y0 = -0.77;
zmin = 4;
zmax = 20;
vg_plot(zg<zmin|zg>zmax) = 10^-9;
i0 = find(zgg==zmin);
iN = find(zgg==zmax);
vg_mask = vg_plot;
count = 0;
% xc = zeros(1,iN-i0+1);
% yc = zeros(1,iN-i0+1);
% zc = zeros(1,iN-i0+1);
% Ic = zeros(1,iN-i0+1);
for iz = i0:iN
    count = count+1;
%     xc(count) = x0;
%     yc(count) = y0;
%     zc(count) = zgg(iz);
    vg1 = squeeze(vg_plot(:,:,iz));
    d = sqrt((xgg(:)-x0).^2+(ygg(:)-y0).^2);
    vg1_search = vg1;
    vg1_search(d>R) = 0;
    vg1_base = max(vg1_search(:));
    vg1_search = vg1_search/vg1_base;
    ww = zeros(size(vg1_search));
    ww(vg1_search>vg_thresh)=1;
    cc = bwconncomp(ww,4);
    object_prop = regionprops(cc);
    if isempty(object_prop)
        vg_mask(:,:,iz) = vg1;
        continue
    end
    object_area = [object_prop.Area];
    [area_max,i_area]=max(object_area);
%    indc = object_prop(i_area).Centroid;
    ind = cc.PixelIdxList;
    indp = ind{i_area};
    mask = zeros(size(vg1));
    mask(indp) = 1;
    vg_mask(:,:,iz) = vg1.*mask;
%     Ic(count) = max(vg1(:).*mask(:));
%     vg_mask_norm(:,:,iz) = vg_mask(:,:,iz)/Ic(count);
%     x0 = interp1(1:size(xgg,2),xgg(1,:),indc(1));
%     y0 = interp1(1:size(ygg,1),ygg(:,1),indc(2));
end
vg_mask(vg_mask==0) = 10^-9; 




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
    v = vg_mask;
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
% xlim([-20 1]);
% ylim([-8 12]);
xlim([-20 1]);
ylim([-12 12]);
zlim([-2 15]);
plot3([0,0],[0,0],[0,4.2],'y','linewidth',4);
hold off;
%az = 76.5;
az = -57+180;
el = 26;
view([az,el])
camlight('headlight');
title(sprintf('%s',covis.sweep.name(7:21)));
h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
xlabel('Easting of COVIS ( m )')
ylabel('Northing of COVIS ( m )')
zlabel('Height above COVIS base ( m )')
end

