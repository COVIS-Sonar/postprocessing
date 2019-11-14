function [grd] = covis_rectgrid_doppler(grd)
% 
% Defines a rectangular grid structure based on covis inputs 
% given by the grd structure.
%
% The bounds are the grid boundaries [xmin xmax ymin ymax zmin zmax].
%
% -----------------
% Version 1.0 - 08/2010 cjones@apl.washington.edu
% Version 1.1 - 05/2011 xupeng_66@hotmail.com
% Version 1.2 - 06/2011 xupeng_66@hotmail.com (added a grid for the radial velocity
% vr)
% Version 1.3 - 09/2018 guangyux@uw.edu (changed to a new grid structure)

if(isfield(grd.bounds,'xmin')); xmin = grd.bounds.xmin; else xmin = 0; end
if(isfield(grd.bounds,'xmax')); xmax = grd.bounds.xmax; else xmax = 0; end
if(isfield(grd.bounds,'ymin')); ymin = grd.bounds.ymin; else ymin = 0; end
if(isfield(grd.bounds,'ymax')); ymax = grd.bounds.ymax; else ymax = 0; end
if(isfield(grd.bounds,'zmin')); zmin = grd.bounds.zmin; else zmin = 0; end
if(isfield(grd.bounds,'zmax')); zmax = grd.bounds.zmax; else zmax = 0; end

if(isfield(grd.spacing,'dx')); dx = grd.spacing.dx; else dx = 0; end
if(isfield(grd.spacing,'dy')); dy = grd.spacing.dy; else dy = 0; end
if(isfield(grd.spacing,'dz')); dz = grd.spacing.dz; else dz = 0; end

if(isfield(grd,'dimensions')) dims = grd.dimensions; else dims = 3; end

x = xmin:dx:xmax; % x coords of the grid
y = ymin:dy:ymax; % y coords of the grid
z = zmin:dz:zmax; % z coords of the grid

% check for empty grid dimensions
if(isempty(x)); x = xmin; end
if(isempty(y)); y = ymin; end
if(isempty(z)); z = zmin; end

% use meshgrid to define the mesh
if(dims == 3) 
    [grd.x, grd.y, grd.z] = meshgrid(x,y,z); % define grid mesh coords
end
if(dims == 2) 
    [grd.x, grd.y] = meshgrid(x,y); % define grid mesh coords
end

% initialize the grid value and weight matrix
if strcmp(grd.type,'doppler velocity')
    grd.vr_cov = zeros(size(grd.x)); % define an empty grid for the covariance-averaged line-of-sight velocity
    grd.vr_vel=zeros(size(grd.x)); % define an empty grid for the velocity-averaged line-of-sight velocity
    grd.std=zeros(size(grd.x)); % define an empty grid for line-of-sight velocity standard deviation
    grd.covar=zeros(size(grd.x)); % define an empty grid for the covariance function
    grd.w = zeros(size(grd.x)); % define an empty grid for weight function
elseif strcmp(grd.type,'intensity')
    grd.I = zeros(size(grd.x)); % define an empty grid for ping-averaged volume backscattering coefficient
    grd.std = zeros(size(grd.x)); % define an empty grid for the standard deviation of volume backscattering coefficient
    grd.w = zeros(size(grd.x)); % define an empty grid for weight function
else
    error('No Doppler grid is found')
end

% set the bounds and size
grd.size = size(grd.x);
grd.axis = [xmin, xmax, ymin, ymax, zmin, zmax];

end

