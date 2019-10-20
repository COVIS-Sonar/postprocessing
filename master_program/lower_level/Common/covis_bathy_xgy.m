% this function is used to reframe the bathymetry of Grotto using the
% bottom of COVIS as the origin.

function [xb,yb,zb,depth_covis,covis_e,covis_n] = covis_bathy_xgy(covis,type,height)

% input: 
% covis: covis structure array
% type: 'diffuse flow', 'imaging', 'doppler'
% height: 1: zb is the height above COVIS ( m )
%         0: zb is the depth ( m )         
% output:
% xb: x-coordinates of grid points ( m )
% yb: y-coordinates of grid points ( m )
% zb: height abvoe COVIS if "height = 1"; depth if "height = 0" ( m )
% depth_covis: depth of COVIS ( m )
% covis_e: UTM East coordinate of COVIS ( m )
% covis_n: UTM North coordinate of COVIS ( m )

if strcmp(type,'imaging')
    grd = covis.grid;
else
    grd = covis.grid{1};
end


% set the offset between the bathymetry data and COVIS data
date_unix = covis.sweep.timestamp{1};
date_mat = unixtime2mat(date_unix);
if date_mat > datenum(2013,6,20,0,0,0);
     coff = 0;
     xoff = 10.02;
     yoff = -0.29;
     zoff = 0;
else
     coff = 0;
     xoff = 8.52;
     yoff = -1.29;
     zoff = 0;
end

% load bathymetry data
file_name = 'EndeavMEF_MBARIAUV1m!.mat'; % MBARIAUV data gridded on a 1-m resolution grid

load(file_name);
if(isfield(covis.sonar.position, 'easting'))
    utm_x0 = covis.sonar.position.easting;
    utm_y0 = covis.sonar.position.northing;
else
    fprintf('No sonar positions available\n');
    return;
end

    
z = Depth;
x = utm_x - (utm_x0); 
y = utm_y - (utm_y0); 

covis_e = utm_x0-xoff;
covis_n = utm_y0-yoff;

% % translate and rotate bathy to fit covis data
R = sqrt(x.^2 + y.^2);
theta = atan2(y, x) + coff;
x = R.*cos(theta) + xoff;
y = R.*sin(theta) + yoff;
z = z + zoff;





%rotate(hb,[0,0,1],theta_off,[0,0,0]);
%rotate(hc,[0,0,1],theta_off,[0,0,0]);

% resample bathy onto a uniform rectangular grid
% this is necessary because utm_x, and utm_y are not necessarily uniform
xmin = grd.bounds.xmin;
xmax = grd.bounds.xmax;
ymin = grd.bounds.ymin;
ymax = grd.bounds.ymax;
dx = grd.spacing.dx;
dy = grd.spacing.dy;

xb = xmin:dx:xmax;
yb = ymin:dy:ymax;
% xb = -40:0.1:10;
% yb = -40:0.1:10;
[xb,yb] = meshgrid(xb,yb);



F = scatteredInterpolant(x(:),y(:),z(:));
zb = F(xb,yb); % gridded height above the bottom of COVIS (m)
depth_covis = F(0,0);
if height == 1
   zb = zb-depth_covis; % height above the bottom of COVIS (m)
end
end

