% This function creates a 2D map of hydrothermal sources on the seafloor 
% within a sector based on a selected acoustic indicator

% version 1.0 by guangyux@uw.edu (Feb 27, 2020)

function covis_diffuse_plot(covis,ind)
% Input:
%  covis: 
%   Matlab structure that includes the gridded data and metadata
%   from lower-level processing of the Diffuse-flow mode data.
%  ind: 
%   acoustic indicator of near-bottom hydrothermal anomalies used for
%   mapping. The options are: ind=1 for 'decorrelation', ind=2 for 
%   'phase variance', and ind=3 for 'amplitude variance'.
% Output:
%  a 2D map of the spatial distribution of hydrothermal sources within a
%  sector based on the selected acoustic indicator. 

% Example:
%  create a map using 'decorrelation' as the acoutic indicator
%  covis_diffuse_plot(covis,1) 

% load gridded data
if ind==1
    grid = covis.grid{2};
    indname = 'Decorr';
    crange = [0.1 0.5];
elseif ind==2
    grid = covis.grid{6};
    indname = 'Phase Variance';
    crange = [0.05 1];
elseif ind==3
    grid = covis.grid{7};
    indname = 'Amplitude Variance';
    crange = [0.05 0.3];
end
xg = grid.x;
yg = grid.y;
v = grid.v;

% load bathymetry data
% swp_name = covis.sweep.name;
% swp_date = datenum(swp_name(7:21),'yyyymmddTHHMMSS');
% if swp_date<=datenum(2019,7,6)
%     bathy_file = sprintf('covis_bathy_2018.mat');
% elseif swp_date<=datenum(2019,11,23)
%     bathy_file = sprintf('covis_bathy_2019a.mat');
% else
%     bathy_file = sprintf('covis_bathy_2019b.mat');
% end
% bathy = load(bathy_file);
% xb = bathy.covis.grid.x;
% yb = bathy.covis.grid.y;
% zb = bathy.covis.grid.v;

% create the map
figure
pcolorjw(xg,yg,v);
axis image;
hold on;
%contour(xb,yb,zb,[-2:0.5:4],'k');
plot(0,0,'.g','markersize',15);
hold off;
xlabel('Easting of COVIS ( m )');
ylabel('Northing of COVIS ( m )');
%title(swp_name);
h = colorbar;
title(h,indname);
ct = cbrewer('seq','Reds',9);
colormap(ct);
caxis(crange);
set(gca,'fontsize',12,'fontname','Times New Roman');


