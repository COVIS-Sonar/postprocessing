function imgfile = covis_imaging_plot(matfile, outputdir, varargin)
%
% Plot covis image grid as isosurfaces
%
% The grid data is loaded from the matfile.
% Plotting parameters are defined with the corresponding
% JSON parameter files (json_file).  If no json file is supplied,
% the default covis_image_plot.json file is used.
% The return string is the figure file name.
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
% ---------------------------
%  Version 1.0 - 10/2010,
%    cjones@apl.washington.edu
%


% Check for other args
p = inputParser;
addParameter(p,'json_file',input_json_path('covis_image_plot_new.json'),@isstring);
parse(p, varargin{:})

json_file = p.Results.json_file;

imgfile = 0;

if(isstruct(matfile))
  covis = matfile;
else
  % pick a mat file, if none given
  if(isempty(matfile))
    error("Matfile %s not specified")
    return
  end

  % check that archive dir exists
  if(~exist(matfile))
      error('Covis .mat file \"%s\" does not exist', matfile);
      return;
  end

  % load the covis gridded data
  load(matfile);
end


if(~isfield(covis,'grid'))
   error('No grid data in covis structure');
   return;
end

% make local copies of the grids
grd = covis.grid; % intensity grid

% parsing the json file
%  which contains all the user supplied parameters
json_str = fileread(json_file);
input = jsondecode(json_str);

if(isfield(input,'verbose'))
    Verbose = input.verbose;
else
    Verbose = 0;
end

if(~isfield(input,'plot_bathy'))
    plot_bathy = 1;
end
plot_bathy = input.plot_bathy;

if(~isfield(input,'isosurface'))
   return;
end
isosurf = input.isosurface

if(Verbose)
    fprintf('Creating isosurface plot for %s\n', covis.sweep.name);
end

if(~isfield(input,'name'))
   input.name = covis.grid.name;
end

if(~isfield(input,'figure'))
   input.figure = 1;
end
fig_num = input.figure;

if(~isfield(input,'visible'))
   input.visible = 1;
end

% creat figure
h = figure(fig_num); clf;
set(gca,'Position',[0.11 0.09 0.775 0.815],'FontSize',30)

% figure visibiliby
if(input.visible == 0)
    set(h, 'Visible', 'off');
end

% grid values
vg = grd.v;
xg = grd.x;
yg = grd.y;
zg = grd.z;

% load bathy
[xb,yb,zb] = covis_bathy(covis, grd);

% plot the bathymetry map
hold on
hb = surf(xb, yb, zb);
shading interp;
colormap('Summer')
contour3(xb, yb, zb, 50,'k');
hold off
%colorbar
%title('COVIS Bathymetry Map for Grotto');
axis equal;

mask = zeros(size(vg));

% mask the data grid using the user inputs
if(isfield(input,'mask'))
    for n=1:length(input.mask)

        m = input.mask(n)

        switch m.type
            case'cylinder'
                x0 = m.x;
                y0 = m.y;
                R = m.radius;
                r = sqrt((xg-x0).^2 + (yg-y0).^2);
                mask(r <= R) = 1;
            case'cone'
                x0 = m.x;
                y0 = m.y;
                z0 = m.z;
                r = sqrt((xg-x0).^2 + (yg-y0).^2);
                R = m.radius;
                H = m.height;
                for m = 1:size(zg,3)
                    zslice = mask(:,:,m);
                    rc = (zg(1,1,m)-z0)*(R/H);
                    zslice(r(:,:,m) <= rc) = 1;
                    mask(:,:,m) = zslice;
                end
            case'parabola'
                x0 = m.x;
                y0 = m.y;
                r0 = m.radius;
                r = sqrt((xg-x0).^2 + (yg-y0).^2);
                mask(r <= r0) = 1;
        end
    end
end

% if mask is empty then unmask everything
if(isempty(find(mask)))
    mask = ones(size(vg));
end

% mask data below the bathy
mask_bathy = 1;
if(mask_bathy)
    dz = 1;
    for m = 1:size(zg,3)
        zslice = ones(size(zg,1),size(zg,2));
        k = find(zg(:,:,m) < (zb+dz));
        zslice(k) = 0;
        mask(:,:,m) = mask(:,:,m) .* zslice;
    end
end

vg(~mask) = nan;

% plot isosurfaces
%figure(h);

hold on;
for n=1:length(isosurf)
   v = vg;
   if(strcmp(isosurf(n).units, 'db'))
      eps = nan;
      m = find(v==0);
      v(m) = eps;  % remove zeros
      v = 10*log10(v);
   end
   surf_value = isosurf(n).value;
   surf_color = isosurf(n).color;
   surf_alpha = isosurf(n).alpha;
   p = patch(isosurface(xg, yg, zg, v, surf_value));
   isonormals(xg, yg, zg, v, p);
   %set(p,'FaceColor','red','EdgeColor','none');
   set(p,'EdgeColor','none','FaceColor',surf_color,'FaceAlpha',surf_alpha);
end
daspect([1 1 1])
hold off;

% set axis
% axis([grd.bounds.xmin grd.bounds.xmax grd.bounds.ymin grd.bounds.ymax ...
%      grd.bounds.zmin grd.bounds.zmax]);
axis([-40 10 -40 10 0 40]);
box on;

% xlabel('Meters','FontSize',20)
% ylabel('Meters','FontSize',20)
zlabel('Distance from COVIS (m)','FontSize',40)

if(isfield(grd,'name') && covis.user.debug)
   str = grd.name;
   str(strfind(str,'_'))='-';
   title(str);
end

% set the view
if(isfield(input,'view'))
   view(input.view.azimuth, input.view.elevation);
else
   view([1,1,1]);
end

%save plot to file
if(~isempty(outputdir))

  % create output dir if it doesn't exist
  if(~exist(outputdir,'dir'))
      warning(['output directory not found, will create one here: ' outputdir])
      mkdir(outputdir);
  end

  % always make a fig file
  figfile = fullfile(outputdir, strcat(input.name, '.fig'));
  if(Verbose) fprintf('Saving figure as %s\n', figfile); end

  if(exist(figfile,'file'))
      fprintf('Warning: overwiting %s\n', figfile);
  end
  saveas(h, figfile);

  % figure type, must be one of formats that the 'saveas' function accepts
  if(isfield(input,'format'))
      type = input.format;
      imgfile = fullfile(outputdir, strcat(input.name,'.',type));
      if(Verbose)
          fprintf('Saving figure %s\n', imgfile);
      end
      if(exist(imgfile,'file'))
          fprintf('Warning: overwiting %s\n', imgfile);
      end
      set(h,'PaperUnits','points','PaperPosition',[0 0 3300 2550],...
          'PaperPositionMode', 'manual','PaperOrientation','portrait','renderer','zbuffer');
      print(h,'-dpng','-r71',imgfile)
  %         saveas(h, imgfile);
  end

end
