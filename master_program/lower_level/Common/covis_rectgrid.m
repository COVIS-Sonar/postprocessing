function [grd] = covis_rectgrid(grd)
% 
% Defines a rectangular grid structure based on covis inputs 
% given by the grd structure.
%
% The bounds are the grid boundaries [xmin xmax ymin ymax zmin zmax].
%
%
% ----------
% This program is free software distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY. You can redistribute it and/or modify it.
% Any modifications of the original software must be distributed in such a 
% manner as to avoid any confusion with the original work.
% 
% Please acknowledge the use of this software in any publications arising  
% from research that uses it.
% 
% -----------------
% Version 1.0 - 08/2010 cjones@apl.washington.edu
%

if(isfield(grd.bounds,'xmin')) xmin = grd.bounds.xmin; else xmin = 0; end
if(isfield(grd.bounds,'xmax')) xmax = grd.bounds.xmax; else xmax = 0; end
if(isfield(grd.bounds,'ymin')) ymin = grd.bounds.ymin; else ymin = 0; end
if(isfield(grd.bounds,'ymax')) ymax = grd.bounds.ymax; else ymax = 0; end
if(isfield(grd.bounds,'zmin')) zmin = grd.bounds.zmin; else zmin = 0; end
if(isfield(grd.bounds,'zmax')) zmax = grd.bounds.zmax; else zmax = 0; end

%xmin = grd.origin.x0 - grd.length/2;
%xmax = grd.origin.x0 + grd.length/2;
%ymin = grd.origin.y0 - grd.width/2;
%ymax = grd.origin.y0 + grd.width/2;
%zmin = grd.origin.z0 - grd.hieght/2;
%zmax = grd.origin.z0 + grd.height/2;
%bounds = [xmin, xmax, ymin, ymax, zmin, zmax];
%grd.bounds = bounds;

if(isfield(grd.spacing,'dx')) dx = grd.spacing.dx; else dx = 0; end
if(isfield(grd.spacing,'dy')) dy = grd.spacing.dy; else dy = 0; end
if(isfield(grd.spacing,'dz')) dz = grd.spacing.dz; else dz = 0; end

if(isfield(grd,'dimensions')) dims = grd.dimensions; else dims = 3; end

x = xmin:dx:xmax; % x coords of the grid
y = ymin:dy:ymax; % y coords of the grid
z = zmin:dz:zmax; % z coords of the grid

% check for empty grid dimensions
if(isempty(x)) x = xmin; end
if(isempty(y)) y = ymin; end
if(isempty(z)) z = zmin; end

% use meshgrid to define the mesh
if(dims == 3) 
    [grd.x, grd.y, grd.z] = meshgrid(x,y,z); % define grid mesh coords
end
if(dims == 2) 
    [grd.x, grd.y] = meshgrid(x,y); % define grid mesh coords
end

% initialize the grid value and weight matrix
grd.v = zeros(size(grd.x)); % define an empty grid value matrix
grd.w = zeros(size(grd.x)); % define an empty grid weight matrix

% set the bounds and size
grd.size = size(grd.v);
grd.axis = [xmin, xmax, ymin, ymax, zmin, zmax];

end

