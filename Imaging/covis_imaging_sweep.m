function matfile = covis_imaging_sweep(swp_file, outputdir, varargin)
%
% Process and grid covis IMAGING sweep data onto a rectangular grid.
%
% The inputs to this function:
%   swp_file - the path to sweep archive file.  Function will handle
%              either a compressed archive (.tar.gz, .zip) or
%              an existing unpacked data directory
%
%   The function will also take additional optional parameters using the
%   input parser method of giving the parameter name, then the value.
%
%    'outputdir' - Directory to store the resulting data structure
%                  as a .mat file
%                  _If not set_ data is not stored to disk
%    'json_file' - Path to a JSON configuration file for the sweep processing
%    'metadata'  - Matlab struct of meta-information to be included in the
%                  struct and .mat file

%   outputdir - directory to save covis structure as mat file
%   json_file - the name of the json input parameter file
%
% The return string is the mat file name that the covis data was saved.
%
% The sweep archive directory contains a set of files for each Covis sweep.
% A sweep is a complete up and down scan if the sonar in elevation.
% The data files in the sweep directory must conform to the file name
% convensions and data formats defined in the document Covis_data.pdf and
% the R7038 data format defined in the document Reson_data_format.html.
% Meta data for the sweep is contained in a sweep.json file.
% The meta data each ping is contianed a the file of the form
% ('rec_7000_%06d.json',ping_number) and the associated binary data file
% name has the form ('rec_7038_%06d.bin',ping_number).
% The binary data files (*.bin) contain element level (quadrature) samples.
% The attitude data for each ping is contained in the 'index.csv' file.
% JSON files are parsed using the function 'parse_json'.
%
% Processing parameters are read from an ascii JSON formatted parameter
% file (json_file) that contains all the information needed to process
% a covis sweep.
% If json_file = [] or zero then a defaults file
% ('input/covis_image.json') is used.
%
% This function returns a structure (covis) that contains all the data and
% meta data for the processed sweep. Processed data is gridded onto an
% evenly spaced rectangular grid stored in the covis.grid structure.
% Gridding is done using the 'l3grid' function and nearest neighbor
% linear interpolation method.
% The data grid is defined in the fixed world coordinate system with the
% sonar in the center of the rectangular coordinate system.
% The y-direction is North, x-direction is East, and z-direction is Up.
% Coordinate transformations from the sensor coordinate system of
% (range, beam azimuth) are done using the function covis_coords().
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
% ----------
%  Version 1.0 - 10/2010,
%    cjones@apl.washington.edu, drj@apl.washington.edu, bemis@rci.rutgers.edu
%  Version 1.1 - 10/2011, cjones@apl.washington.edu
%

global Verbose;

% Check for other args
p = inputParser;
addParameter(p,'json_file',input_json_path('covis_image.json'),@isstring);
addParameter(p,'metadata',0,@isstruct);
parse(p, varargin{:})

% Extract a COVIS archive, if it hasn't been unpacked already
[swp_path, swp_name] = covis_extract(swp_file, '');
swp_dir = fullfile(swp_path, swp_name);

%% On error, return matfile = ''
matfile = '';
covis = struct;

% Create MAT output filename; check if it exists
if(~isempty(outputdir))
    matfile = char(fullfile(outputdir, strcat(swp_name, '.mat')));
    if exist(matfile,'file')
      fprintf('Warning: not overwiting %s\n', matfile);
      return ;
  end
end

fprintf("a\n")

% parse sweep.json file in data archive
swp_file = fullfile(swp_dir, 'sweep.json');
if(~exist(swp_file))
    fprintf('sweep.json file does not exist at %s\n', swp_file);
    return;
end

json_str = fileread(swp_file);
swp = jsondecode(json_str);

% set sweep path and name
swp.path = swp_path;
swp.name = swp_name;

fprintf("b\n")


% parsing the json input file for the user supplied parameters
json_file = p.Results.json_file
fprintf('Using sweep config file %s\n', json_file)

% check that json input file exists
if ~exist(json_file,'file')
    fprintf('JSON input file %s does not exist\n', json_file);
    return;
end
json_str = fileread(json_file);
covis = jsondecode(json_str);
if(strcmpi(covis.type, 'imaging') == 0)
    fprintf('Incorrect covis input file type\n');
    return;
end


Verbose = covis.user.verbose;
Debug_Plot = covis.user.debug;

% define a 3D rectangular data grid
[covis.grid] = covis_rectgrid(covis.grid);
covis.grid.name = swp.name;

% set local copies of covis structs
grd = covis.grid;
usr = covis.user;
pos = covis.sonar.position;
dsp = covis.processing;
bfm = covis.processing.beamformer;
cal = covis.processing.calibrate;
filt = covis.processing.filter;

fprintf("c\n")


% check outpath
if(~isfield(usr,'outpath'))
    usr.outpath = 'output';  % default
    if(Verbose) fprintf('Setting user output path: %s\n',usr.outpath); end;
end

% check the grid
% Set the type of grid value (intensity or complex).
if(~isfield(grd,'type'))
    grd.type = 'intensity';  % default grid type
    if(Verbose) fprintf('Setting grid type: %s\n',grd.type); end;
end

% Set the type of beamforming (fast, fft, ...)
if(~isfield(bfm,'type'))
    bfm.type = 'fast';
    if(Verbose) fprintf('Setting beamform type: %s\n',bfm.type); end;
end
bfm_type = bfm.type;

% Set compass declination
if(~isfield(pos,'declination'))
    pos.declination = 18.0;
end

% calibration parameters
if(~isfield(cal,'mode'))
    cal.mode = 'VSS'; % 'VSS' or 'TS'
end

% set position of sonar
% should use covis.position info, but for now ...
if(~isfield(pos,'altitude'))
    pos.altitude = 4.2;
end
alt = pos.altitude;
origin = [0,0,alt]; % sonar is always at (0,0,0) in world coords

% directory list of *.bin file
file = dir(fullfile(swp_dir, '*.bin'));
nfiles = length(file);

% Read index file
% The index file contains the sweep parameters:
% ping,seconds,microseconds,pitch,roll,yaw,kPAngle,kRAngle,kHeading
% in csv format
ind_file = 'index.csv';
% should check if file exists
if(Verbose)
    fprintf('Parsing %s\n', fullfile(swp_dir, ind_file));
end;
csv = csvread(fullfile(swp_dir, ind_file), 1, 0);
if(size(csv,1) ~= nfiles)
    fprintf('index size and file number mismatch');
end
if(Verbose > 2)
    fprintf('finished parsing: %s\n',datestr(now))
end

% save index data in ping structure
for n=1:nfiles
    png(n).num = csv(n,1);
    png(n).sec = (csv(n,2) + csv(n,3)/1e6);
    % save rotator angles
    png(n).rot.pitch = csv(n,4);
    png(n).rot.roll = csv(n,5);
    png(n).rot.yaw = csv(n,6)';
    % save tcm6 angles
    png(n).tcm.kPAngle = csv(n,7)';
    png(n).tcm.kRAngle = csv(n,8)';
    png(n).tcm.kHeading = (pos.declination + csv(n,9))';
end

% find the number of bursts in the sweep and number of pings in each burst
% by checking when the elevation changes
[ burst ] = covis_parse_bursts( png );

% range of elev steps to process
elev_start = (covis.processing.bounds.pitch.start);
elev_stop = (covis.processing.bounds.pitch.stop);
range_start = (covis.processing.bounds.range.start);
range_stop = (covis.processing.bounds.range.stop);

% loop over bursts
nbursts = length(burst);
for nb = 1:nbursts

    % check elevation
    if((burst(nb).elev < elev_start) | (burst(nb).elev > elev_stop))
        nb_start = nb + 1;
        continue;
    end

    if(Verbose)
        fprintf('Burst %d: Elevation %0.2f\n', nb, burst(nb).elev);
    end;

    npings = burst(nb).npings;

    % check that there's enough pings in burst
    if((npings < 2) & strfind(dsp.ping_combination.mode,'diff'))
        fprintf('Not enough pings in burst\n');
        continue;
    end

    % loop over pings in a burst, read data and hold onto it
    for np = 1:npings

        n = burst(nb).start_ping + (np-1); % ping number

        % parse filenames
        bin_file = file(n).name;
        [type,ping_num] = strread(bin_file,'rec_%d_%d.bin');
        if(ping_num ~= n)
            fprintf('Ping number mismatch: %d\n', n);
        end
        json_file = sprintf('rec_7000_%06d.json',ping_num);

        % read ping meta data from fson file
        json_str = fileread(fullfile(swp_dir, json_file));
        json = jsondecode(json_str);
        png(n).hdr = json.hdr;

        % define sonar attitude (use TCM6 angles)
        elev = (pi/180) * png(n).tcm.kPAngle;
        roll = (pi/180) * png(n).tcm.kRAngle;
        yaw = (pi/180) * png(n).tcm.kHeading;

        if(Verbose > 1)
            fprintf(' %s: elev %0.2f, roll %0.2f, yaw %0.2f\n', bin_file, ...
                elev*180/pi, roll*180/pi, yaw*180/pi);
        end;

        % read raw element quadrature data
        [hdr, data] = covis_read(fullfile(swp_dir, bin_file));

        if(np == 1)
            monitor = data;
        end

        % Correct phase using first ping as reference
        data = covis_phase_correct(png(n), monitor, data);

        % Apply Filter to data
        [data, filt, png(n)] = covis_filter(data, filt, png(n));

        % define beamformer parameters
        bfm.fc = png(n).hdr.xmit_freq;
        bfm.c = png(n).hdr.sound_speed;
        bfm.fs = png(n).hdr.sample_rate;
        bfm.first_samp = hdr.first_samp + 1;
        bfm.last_samp = hdr.last_samp;
        % beamform
        [bfm, bf_sig(:,:,np)] = covis_beamform(bfm, data);

        % save ping in the burst array
        bf_sig(:,:,np) = covis_calibration(bf_sig(:,:,np), bfm, png(n), cal);

        % plot each ping
        %range = bfm.range;
        %azim = bfm.angle;
        %xs = range*sin(azim);
        %ys = range*cos(azim);
        %figure(1);
        %subplot(3,2,np);
        %pcolor(xs,ys,20*log10(abs(bf_sig(:,:,np)))); shading flat; axis equal;
        %caxis([-90 -30]);
        %axis([5 25 25 45]);
        %%colorbar;
        %pause;

    end

    % calc the value to grid
    if(strfind(dsp.ping_combination.mode,'diff'))
        % mean of the abs^2 of diff between pings
        v = mean(abs(0.5*diff(bf_sig,1,3)).^2, 3);
    else
        % mean of the intensity of the pings
        v = mean(abs(bf_sig).^2, 3);
    end

    n = burst(nb).start_ping; % ping number
    % define sonar attitude (use TCM6 angles)
    elev = (pi/180) * png(n).tcm.kPAngle;
    roll = (pi/180) * png(n).tcm.kRAngle;
    yaw = (pi/180) * png(n).tcm.kHeading;

    % transform sonar coords into world coords
    range = bfm.range;
    azim = bfm.angle;
    [xv, yv, zv] = covis_coords(origin, range, azim, yaw, roll, elev);

    % grid the ping data
    [grd.v,grd.w] = l3grid(xv,yv,zv,v,grd.x,grd.y,grd.z,grd.v,grd.w);

    if(Verbose > 2)
        fprintf('finished gridding: %s\n',datestr(now))
    end

    if(Debug_Plot)
        figure(1); clf;
        %pcolor(xv,yv,v); shading flat; axis equal; caxis([0 1]);
        surf(xv,yv,zv,v); shading flat; axis equal;
        %surf(xv,yv,zv,10*log10(abs(bf_sig).^2)); shading flat; axis equal;
        %axis([-50 50 -50 50 -10 40]);
        if(isfield(usr,'view'))
            view(usr.view.azimuth, usr.view.elevation);
        else
            view([1,1,1]);
        end
        str = sprintf('Ping Number %d: Elevation %.2f, Roll %.2f, Yaw %.2f\n', ...
            png(n).hdr.ping_num, elev*180/pi, roll*180/pi, yaw*180/pi);
        title(str);
        colorbar
        refresh;
        %pause;
    end;

end   % End loop over bursts

% normalize the grid with the grid weights
n = find(grd.w);
grd.v(n) = grd.v(n)./grd.w(n);

% save local copies of covis structs
covis.sweep = swp;
covis.grid = grd;
covis.ping = png;
covis.sonar.position = pos;
covis.processing.beamformer = bfm;
covis.processing.calibrate = cal;
covis.processing.filter = filt;
covis.burst = burst;

fprintf("m\n")


% save covis structure in a mat file for later use
if(~isempty(matfile))
  fprintf("Saving results to %s\n", matfile)

    if(~exist(outputdir,'dir'))
        mkdir(outputdir);
    end

    covis.metadata = p.Results.metadata

    save(matfile,'covis');
end

fprintf("n\n")


end
