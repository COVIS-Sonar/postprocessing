function imgfile = covis_diffuse_plot(matfile, outputdir, varargin)
%
%function [h, covis, imgfile] = covis_diffuse_plot(swp_path, covis, varargin)
%
% Plot covis diffuse grid
%
% Inputs:
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
%    cjones@apl.washington.edu, drj@apl.washington.edu

% normalize the data for now
%maxv = max(grd.v(:));
%v = grd.v / maxv;

% swp_path is added as a new input parameter, which and outpath in json file will be used provide
% directory for output files. Previously only outpath in json file is used to
% determine directory for output files.--added by yingsong in Oct 2011.

% do_plot_bathy = 1;

% Check for other args
p = inputParser;
addParameter(p,'json_file',input_json_path('covis_diffuse_plot.json'),@isstring);
parse(p, varargin{:})

json_file = p.Results.json_file;

if(isstruct(matfile))
  covis=matfile
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
  load(matfile)

end


json_str = fileread(json_file);
input = parse_json(json_str);

% make local copies of the grids
for n=1:length(covis.grid)
    switch lower(covis.grid(n).type)
        case {'decorrelation'}
            grd = covis.grid(n); % decorr grid
        case 'intensity'
            I_grd = covis.grid(n); % intensity grid
        case 'decorrelation intensity'
            dI_grd = covis.grid(n); % decorr intensity grid
        otherwise disp('Unknown grid type.')
    end
end

if(~isfield(covis.user,'figure'))
    covis.user.figure = 1;
end
fig_num = covis.user.figure;

% load bathy grid
[xb,yb,zb] = covis_bathy(covis, grd);

vb = nan*ones(size(zb));

% drape corr grid over bathy

x = grd.x;
y = grd.y;

% grid bathy onto uniform corr grid in x and y
z = griddata(xb, yb, zb, x, y);

R3 = sqrt(x.^2 + y.^2 + z.^2); % range to 3D point in space
R2 = sqrt(x.^2 + y.^2); % range to flat bottom
theta = atan2(y, x);
v2 = grd.v;
v2(v2==0) = nan;

%v = griddata(R2, theta, v2, R3, theta, 'nearest');
v = griddata(R2, theta, v2, R3, theta, 'linear');

% w = zeros(size(x));
% v = zeros(size(x));
% [v, w] = l2grid(R2, theta, grd.v, R3, theta, v, w);
% % normalize the grid with the grid weights
% n = find(w);
% v(n) = v(n)./w(n);

v(v==0) = nan;

% plot the bathymetry map
h = figure(fig_num); clf;
set(gca,'Position',[0.05 0.1 0.8 0.85],'FontSize',24)
hold on
hb = surf(x, y, z, v);
shading flat;
colormap('hot')
contour3(x, y, z, 0:1:20, 'w');
hold off
%colorbar
%title('COVIS Bathymetry Map for Grotto');
axis equal;

if strcmp(grd.type, 'decorrelation')
    caxis([0 0.5]);
    %axis([-50 50 -50 50 -10 40]);
    titleText = 'Average Decorrelation';
elseif strcmp(grd.type, 'decorrelation intensity')
    caxis([-30 0]);
    %axis([-50 50 -50 50 -10 40]);
    titleText = 'Decorrelation Intensity (dB re max)';
else
    caxis([-30 0]);
    %axis([-50 50 -50 50 -10 40]);
    titleText = 'Intensity (dB re max)';
end

%if(isfield(grd,'name'))
%   str = grd.name;
%   str(strfind(str,'_'))='-';
%   title(str);
%end

% zlabel('Meters','FontSize',40);
xlabel('Distance from COVIS (m)','FontSize',40);

if(isfield(covis.user,'view'))
   view(input.view.azimuth, input.view.elevation);
else
    view([0,0,1]);
end

cb_axes = axes('position',[0.9 0.1 0.05 0.75],'FontSize',40);
colorbar(cb_axes)
title(cb_axes,titleText,'FontSize',40);

% figure type, must be one of formats that the 'saveas' function accepts
if(~isfield(input,'format'))
   input.format = 'fig';
end
type = input.format;

%save plot to file
if(~isempty(outputdir))
    % swp_path is also applied to determine directory for output files
    %input.outpath=strcat(swp_path,input.outpath);

   if ~exist(outputdir,'dir')
       if ~isempty(outputdir)
           mkdir(outputdir);
       end
   end

   imgfile = fullfile(outputdir, strcat(grd.name,'.',type));

   if(exist(imgfile,'file'))
      fprintf('Warning: overwiting %s\n', imgfile);
   end

   set(h,'PaperUnits','points','PaperPosition',[0 0 3300 2550],...
        'PaperPositionMode', 'manual','PaperOrientation','portrait','renderer','zbuffer');
   print(h,'-dpng','-r71',imgfile)
   %saveas(h, imgfile);

end

if covis.user.debug
    figure(fig_num+1); clf;
    pcolor(x, y, v2);
    shading flat; axis equal;

    if strcmp(grd.type, 'decorrelation')
        caxis([0 1]);
        %axis([-50 50 -50 50 -10 40]);
        title('Average Decorrelation')
    elseif strcmp(grd.type, 'decorrelation intensity')
        caxis([-30 0]);
        %axis([-50 50 -50 50 -10 40]);
        title('Decorrelation Intensity (dB re max)')
    else
        caxis([-30 0]);
        %axis([-50 50 -50 50 -10 40]);
        title('Intensity (dB re max)')
    end
end

end
