function matfile = covis_diffuse_sweep(swp_file, outputdir, varargin)
%
% Process and grid covis DIFFUSE sweep data onto a rectangular grid.
%
% The inputs to this function:
%  swp_path - the path to the sweep directory
%  swp_name - the name of the sweep directory
%  json_file - the name of the json parameter file
%
% Processing parameters are read from an ascii JSON formatted parameter
% file (json_file) that contains all the information needed to process
% a covis sweep.
%
% The sweep directory (swp_dir = fullfile(swp_path, swp_name)) contains
% a set of files for each Covis sweep. A diffuse sweep is a complete
% set of pings at a fixed sonar eleveation.
% If json_file = [] or zero then a defaults file is used,
% named ('input/covis_diffuse.json').
%
% The sweep directory contains a set of data files for each ping.
% The files in the sweep directory must conform to the file name
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
% This function returns a structure (covis) that contains all the data and
% meta data for the processed sweep. Processed data is gridded onto an
% evenly spaced rectangular grid stored in the covis.grid structure.
%
% Gridding is done using the 'l2grid' function and nearest neighbor
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

% swp_path is also applied to determine directory for output files--added by yingsong in Oct 2011

global Verbose;

% Check for other args
p = inputParser;
addParameter(p,'json_file',input_json_path('covis_diffuse.json'),@isstring);
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
      return
  end
end

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

% parsing the json input file for the user supplied parameters
json_file = p.Results.json_file
fprintf('Using sweep config file %s\n', json_file)

% check that json input file exists
if ~exist(json_file,'file')
    fprintf('JSON sweep config file %s does not exist\n', json_file);
    return;
end
json_str = fileread(json_file);
covis = jsondecode(json_str);
if(strcmpi(covis.type, 'diffuse') == 0)
    fprintf('JSON sweep config file of incorrect type \"%s\"\n', covis.type);
    return;
end

Verbose = covis.user.verbose;
Debug_Plot = covis.user.debug;

disp(covis.grid)

% define a 2D rectangular data grid
newgrid = struct(covis_rectgrid(covis.grid(1))) % corr grid
newgrid.name = swp.name
for n=2:length(covis.grid)
  grd = covis_rectgrid(covis.grid(n)) % corr grid
  grd.name = swp.name;
  newgrid(n) = grd
end

covis.grid = newgrid

% set local copies of covis structs
usr = covis.user;
pos = covis.sonar.position;
dsp = covis.processing;
bfm = covis.processing.beamformer;
cal = covis.processing.calibrate;
filt = covis.processing.filter;
cor = covis.processing.correlation;

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
    cal.mode = 'TS'; % 'VSS' or 'TS'
end
if(~isfield(cal,'filt'))
    cal.filt = 1;         % 1 for low-pass filtering, 0 for no filtering
    cal.filt_bw = 2.0;     % Bandwdith is filt_bw/tau
end

% sonar height off bottom
if(~isfield(pos,'altitude'))
    pos.altitude = 4.2;
end
height = pos.altitude;  % sonar height off bottom

% set position of sonar
% should use covis.position info, but for now ...
origin = [0,0,height];

% correlation range window size [number of samples]
if(~isfield(cor,'window_size'))
    cor.window_size = 0.001;
end
window_size = cor.window_size;

% correlation window overlap [number of samples]
if(~isfield(cor,'window_overlap'))
    cor.window_overlap = 0.4;
end
window_overlap = cor.window_overlap;

% Threshold for windowing image re max intensity
windthresh = cor.windthresh;

% correlation lag number of samples
nlag = cor.nlag;

% directory list of *.bin file
file = dir(fullfile(swp_dir, '*.bin'));
nfiles = length(file);

if nfiles <= nlag
    warning('COVIS:covis_diffuse_sweep:badTar',['Warning: covis tar file ' swp_dir ' is missing data files. Processing of this data set cannot be completed at this time.']);
    covis = struct([]);
    matfile = '';
    return;
end

% read index file
% The index file contains the sweep parameters:
% ping,seconds,microseconds,pitch,roll,yaw,kPAngle,kRAngle,kHeading
% in csv format
ind_file = 'index.csv';
try
    csv = csvread(fullfile(swp_dir, ind_file), 1, 0);
catch
    warning('COVIS:covis_diffuse_sweep:badIndexFile',['Warning: in covis tar file ' swp_dir ' the index file cannot be read. Processing of this data set cannot be completed at this time.']);
    covis = struct([]);
    matfile = '';
    return;
end
if(size(csv,1) ~= nfiles)
    warning('COVIS:covis_diffuse_sweep:badIndexFile',['Warning: in covis tar file ' swp_dir ' the index size and file number do not match. Processing of this data set cannot be completed at this time.']);
    covis = struct([]);
    matfile = '';
    return;
end

if(Verbose)
    fprintf('Parsing %s\n', fullfile(swp_dir, ind_file));
end;

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


% loop over ping files
for n=1:nfiles

    % parse filenames
    bin_file = file(n).name;
    [type,ping_num] = strread(bin_file,'rec_%d_%d.bin');
    json_file = sprintf('rec_7000_%06d.json',ping_num);

    %if(Verbose)
    %    fprintf('Parsing %s\n', fullfile(swp_dir, json_file));
    %end;

    % read ping meta data from fson file
    json_str = fileread(fullfile(swp_dir, json_file));
    json = parse_json(json_str);
    png(n).hdr = json.hdr;

    if (n == 1)
        if(Verbose > 1)
            png(1).hdr  % View essential parameters
        end
        fsamp = png(1).hdr.sample_rate;
        cwsize = round(fsamp*window_size);
        cwovlap = round(window_overlap*cwsize);
    end

    pitch = (pi/180) * png(n).tcm.kPAngle;
    roll = (pi/180) * png(n).tcm.kRAngle;
    yaw = (pi/180) * png(n).tcm.kHeading;

    if(Verbose > 1)
        fprintf('Reading %s: pitch %f, roll %f, yaw %f\n', bin_file, ...
            pitch*180/pi, roll*180/pi, yaw*180/pi);
    end;

    % read raw element quadrature data
    [hdr, data(:,:,n)] = covis_read(fullfile(swp_dir, bin_file));

end   % End loop on pings


% Loop again on pings to form ping-ping correlation
for n = 1:nfiles-nlag

    % Take two pings separated by lag
    data1 = data(:,:,n);
    data2 = data(:,:,n+nlag);

    % Correct phase
    data2 = covis_phase_correct(png(n), data1, data2);

    % define beamformer parameters
    bfm.fc = png(n).hdr.xmit_freq;
    bfm.c = png(n).hdr.sound_speed;
    bfm.fs = png(n).hdr.sample_rate;
    bfm.first_samp = hdr.first_samp + 1;
    bfm.last_samp = hdr.last_samp + 1;

    % beamform the quadrature data
    [bfm, bf_sig1] = covis_beamform(bfm, data1);
    [bfm, bf_sig2] = covis_beamform(bfm, data2);

    % apply calibration to beamformed data
    bf_sig1 = covis_calibration(bf_sig1, bfm, png(n), cal);
    bf_sig2 = covis_calibration(bf_sig2, bfm, png(n), cal);

    % Apply Filter to data
    [bf_sig1, filt, png(n)] = covis_filter(bf_sig1, filt, png(n));
    [bf_sig2, filt, png(n)] = covis_filter(bf_sig2, filt, png(n));

    % define sonar attitude
    pitch = (pi/180) * png(n).tcm.kPAngle;
    roll = (pi/180) * png(n).tcm.kRAngle;
    yaw = (pi/180) * png(n).tcm.kHeading;

    % transform sonar coords into world coords
    range = bfm.range;
    azim = bfm.angle;

    % Correlate pings
    %  rc is the range of the center of the corr bin
    [cor, I, rc] = covis_cor_win(bf_sig1, bf_sig2, range, cwsize, cwovlap);

    % window in range
    ind = find(rc > height);
    rc = rc(ind);
    I = I(ind,:);
    cor = cor(ind,:);

    % save the quanities of interest
    cor_av(:,:,n) = abs(cor);
    I_av(:,:,n) = abs(I);
    decor_I_av(:,:,n) = I.*(1.0 - abs(cor));

    clear data1 data2

end  % End loop on ping number, n

% average
cor_av = mean(cor_av,3);
I_av = mean(I_av,3);
decor_I_av = mean(decor_I_av,3);

% Find location of range bins on a flat bottom
r = sqrt(height^2 - rc.^2); % range along bottom
for m = 1:length(azim)
    xv(:,m) = rc*cos(-azim(m) + yaw);
    yv(:,m) = rc*sin(-azim(m) + yaw);
    zv(:,m) = zeros(length(rc),1);
end

% grid the data, loop over different grid types
for n = 1:length(covis.grid)

    grd = covis.grid(n);

    switch lower(grd.type)
        case 'decorrelation'
            % Window decorrelation based on intensity
            Inorm = I_av/max(max(I_av));
            window = 1 - exp(-(Inorm/windthresh+0.0001).^4);
            v = window.*(1.0 - cor_av);
            %v = (1.0 - cor_av);
        case 'intensity'
            % Normalize to max of unity
            v = I_av/max(max(I_av));
        case 'decorrelation intensity'
            % Normalize to max of unity
            v = decor_I_av/max(max(decor_I_av));
        otherwise disp('Unknown grid type.')
    end

    % grid the value
    [grd.v, grd.w] = l2grid(xv, yv, v, grd.x, grd.y, grd.v, grd.w);

    % normalize the grid with the grid weights
    m = find(grd.w);
    grd.v(m) = grd.v(m)./grd.w(m);

    covis.grid(n) = grd;

end

% save local copies of covis structs
covis.sweep = swp;
covis.beamformer = bfm;
covis.ping = png;
covis.calibrate = cal;

% save grd for later use
if(~isfield(usr,'outpath'))
   % swp_path is also applied to determine directory for output files
    usr.outpath = strcat(swp_path,'/output');  % default
   % usr.outpath = '/output';  % default
    if(Verbose) fprintf('Setting user output path: %s\n',usr.outpath); end;
else
    usr.outpath = strcat(swp_path,usr.outpath);
    % usr.outpath = usr.outpath;
end

if(~exist(usr.outpath,'dir'))
    if ~isempty(usr.outpath)
        mkdir(usr.outpath);
    end
end

if(~isempty(matfile))
  fprintf("Saving results to %s\n", matfile)

    if(~exist(outputdir,'dir'))
        mkdir(outputdir);
    end

    covis.metadata = p.Results.metadata

    save(matfile,'covis');
end

if(Debug_Plot)

    for n = 1:length(covis.grid)
        grd = covis.grid{n};
        figure; clf;
        pcolor(grd.x,grd.y,grd.v); shading flat; axis equal;
        switch lower(grd.type)
            case 'decorrelation'
                %caxis([0 1]);
                %axis([-50 50 -50 50 -10 40]);
                str = sprintf('Average Decorrelation Intensity of Pings 1-%d: Lag %.2f\n', ...
                    png(n).hdr.ping_num, png(n).hdr.ping_period);
                %        title(str);
                title('Average Decorrelation')
            case 'decorrelation intensity'
                %caxis([-30 0]);
                %axis([-50 50 -50 50 -10 40]);
                str = sprintf('Average Decorrelation Intensity of Pings 1-%d: Lag %.2f, Pitch %.2f, Roll %.2f, Yaw %.2f\n', ...
                    png(n).hdr.ping_num, png(n).hdr.ping_period, pitch*180/pi, roll*180/pi, yaw*180/pi);
                %        title(str);
                title('Decorrelation Intensity (dB re max)')
            case 'intensity'
                %axis([-50 50 -50 50 -10 40]);
                str = sprintf('Average Intensity of Pings 1-%d: Lag %.2f, Pitch %.2f, Roll %.2f, Yaw %.2f\n', ...
                    png(n).hdr.ping_num, png(n).hdr.ping_period, pitch*180/pi, roll*180/pi, yaw*180/pi);
                %        title(str);
                title('Intensity (dB re max)')

        end
        colorbar
        refresh;

    end
end

end
