% This program is used to combine the gridded bathymetry data recorded in
% all three sectors of COVIS's field-of-view and plot the resulting seafloor
% topography

% version 1.0 by guangyux@uw.edu (Oct 18, 2019)

% load gridded bathymetry data
data_dir = 'F:\COVIS\Axial\COVIS_data\processed\Bathymetry\2019';
data_name1 = 'COVIS-20190706T040730-bathy1'; % Sector 1 gridded bathymetry
data_name2 = 'COVIS-20190706T042134-bathy2'; % Sector 2 gridded bathymetry
data_name3 = 'COVIS-20190706T043548-bathy3'; % Sector 3 gridded bathymetry
bathy1 = load(fullfile(data_dir,data_name1));
bathy2 = load(fullfile(data_dir,data_name2));
bathy3 = load(fullfile(data_dir,data_name3));

% combine gridded bathymetry
xg = bathy1.covis.grid.x;
yg = bathy1.covis.grid.y;
zb(:,:,1) = bathy1.covis.grid.v;
w(:,:,1) = bathy1.covis.grid.w;
zb(:,:,2) = bathy2.covis.grid.v;
w(:,:,2) = bathy2.covis.grid.w;
zb(:,:,3) = bathy3.covis.grid.v;
w(:,:,3) = bathy3.covis.grid.w;
zb_sum = nansum(zb.*w,3);
w_sum = nansum(w,3);
zb_out = zb_sum;
zb_out(w_sum~=0) = zb_sum(w_sum~=0)./w_sum(w_sum~=0);
zb_out(zb_out==0) = nan;
covis.grid.x = xg;
covis.grid.y = yg;
covis.grid.v = zb_out;
covis.grid.w = w_sum;

% save combined gridded bathymetry
out_dir = data_dir;
out_name = [data_name1(1:14),'_bathy_combine'];
save(fullfile(out_dir,out_name),'covis');

% plot combined gridded bathymetry
figure
pcolorjw(xg,yg,zb_out);
axis image;
colormap('jet');
caxis([-1 5]);


