% This function processes the raw data recorded in COVIS's
% bathymetry mode to generate gridded data that is saved in the 'covis'
% structure array

% version 1.0 by guangyux@uw.edu (Oct 19, 2019)
%  --based on the original code written by Chris Jones in 2010
% verison 2.0 by guangyux@uw.edu (Feb 27, 2020)
%  --more adjustable parameters are added to the metadata in the output
%    structure
%  --version number of the code is added to the metadata in the output structure

function covis = covis_bathy_sweep(swp_path, swp_name, json_file, fig)
% Input:
% swp_path: raw data directory
% swp_name: name of raw data sweep
% json_file: name of the json file that includes the metadata needed for
% processing. If 0, use the default json file
% fig: set fig to 1 to plot gridded data. Plotting is muted otherwise.
% Output:
% covis: this is the Matlab structure array that includes the gridded data
% and metadata

% Example:
% swp_path = 'C:\COVIS\Axial\COVIS_data\raw\Bathy\2019';
% swp_name = 'COVIS-20191227T155911-bathy1';
% json_file = 0;
% fig = 1;

%% Initialization
% define cali_year for selection between 2010 and 2018 calibration files

% sonar's central yaw and heading
year = swp_name(7:10);
central_yaw = 135; % central yaw motor reading
if strcmp(year,'2018')
    central_head = 289; % sonar's central magnetic heading measured in 2018 ( degree )
else
    central_head = 279; % sonar's central magnatic heading measured in 2019 ( degree )
end

% calibration parameters
T = 2.41; % temperataure (degree C)
S = 34.53; % salinity (PSU)
pH = 7; % pH
lat = 46; % latitude ( degree N )
depth =1544; % depth ( m )
p = gsw_p_from_z(-depth,lat);
c = gsw_sound_speed_t_exact(S,T,p); % sound speed ( m/s )


% signal-to-noise (SNR) parameters
noise_floor = 0.1970; % rms noise floor
snr_thresh = 60; % SNR threshold (dB)

% OSCFAR parameters
clutterp = 65; % start with the 65% quantile as the clutter threshold
scrthreshold = 20; % dB threshold for signal-to-clutter ratio


% START OF MAIN PROGRAM
swp_dir = fullfile(swp_path, swp_name);

% check that archive dir exists
if(~exist(swp_dir,'dir'))
    error('Sweep directory %s does not exist\n',swp_dir);
end

% parse sweep.json file in data archive
swp_file = 'sweep.json';
if(~exist(fullfile(swp_dir, swp_file),'file'))
    error('sweep.json file does not exist in %s',swp_dir);
end
json_str = fileread(fullfile(swp_dir, swp_file));
swp = parse_json(json_str);

% set sweep path and name
swp.path = swp_path;
swp.name = swp_name;

% parsing the json input file for the user supplied parameters
if isempty(json_file) || all(json_file == 0)
    % default json input file
    json_file = find_input_file('covis_bathy.json');
end
% check that json input file exists
if ~exist(json_file,'file')
    error('JSON input file does not exist');
end

% parse json file
json_str = fileread(json_file);
covis = parse_json(json_str);
if(strcmpi(covis.type, 'bathy') == 0)
    error('Incorrect covis input file type. Looking for: bathy, Current type: %s\n',covis.type);
end

% Define global variable Verbose and read its value from the json file
global Verbose;
Verbose = covis.user.verbose;


% define a 2D rectangular data grid
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


% check the grid
% Set the type of grid Ialue (intensity or complex).
if(~isfield(grd,'type'))
    grd.type = 'intensity';  % default grid type
    if(Verbose)
        fprintf('Setting grid type: %s\n',grd.type);
    end
end

% Set the type of beamforming (fast, fft, ...)
if(~isfield(bfm,'type'))
    bfm.type = 'fast';
    if(Verbose)
        fprintf('Setting beamform type: %s\n',bfm.type);
    end
end


% Set compass declination
if(~isfield(pos,'declination'))
    pos.declination = 16.0;
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

% check if there are missing json files
json_file = dir(fullfile(swp_dir, '*.json'));
if length(json_file)<nfiles
    fprintf('Warning: there are %d json files missing in the sweep', nfiles-length(json_file));
elseif length(json_file)>nfiles
    fprintf('Warning: there are %d bin files missing in the sweep', length(json_file)-nfiles);
end

% read ping meta data from json file
json_str = fileread(fullfile(swp_dir, json_file(1).name));
json = parse_json(json_str);
ind = 1;
while abs(json.hdr.sample_rate-34483)>1 && ind<=length(json_file)
    ind = ind+1;
    json_str = fileread(fullfile(swp_dir, json_file(ind).name));
    json = parse_json(json_str);
end

% Read index file
% The index file contains the sweep parameters:
% ping,seconds,microseconds,pitch,roll,yaw,kPAngle,kRAngle,kHeading
% in csv format
ind_file = 'index.csv';
% should check if file exists
if(Verbose)
    fprintf('Parsing %s\n', fullfile(swp_dir, ind_file));
end
csv = csvread(fullfile(swp_dir, ind_file), 1, 0);
csv = csv(csv(:,1)~=0,:);


% save index data in ping structure
for n=1:size(csv,1)
    png(n).num = csv(n,1);
    png(n).sec = (csv(n,2) + csv(n,3)/1e6);
    png(n).rot_pitch = csv(n,4);
    png(n).sen_pitch = csv(n,7);
    png(n).rot_roll = csv(n,5)/6;
    png(n).sen_roll = csv(n,8);
    if csv(n,6) == 0
        png(n).rot_yaw = 0;
    else
        png(n).rot_yaw = csv(n,6)-central_yaw;
    end
    png(n).sen_head = csv(n,9);
    png(n).hdr = json.hdr;
end

% group the pings transmitted at the same elevation angles into bursts
[ burst ] = covis_parse_bursts_pitch( png );


% range of elev angles to process
elev_start = (covis.processing.bounds.pitch.start);
elev_stop = (covis.processing.bounds.pitch.stop);

%% Main program
% loop over bursts
bad_ping = zeros(0);
bad_ping_count = 0;
nbursts = length(burst);
burst_count = 0;
for nb = 1:nbursts

    % check elevation
    if((burst(nb).pitch < elev_start) || (burst(nb).pitch > elev_stop))
        continue;
    end

    if(Verbose)
        fprintf('Burst %d: pitch %0.2f\n', nb, burst(nb).pitch);
    end

    npings = burst(nb).npings;

    % check that there's enough pings in burst
    if((npings < 2))
        fprintf('Warning: not enough pings in burst\n');
        continue;
    end

    % loop over pings in a burst, read data and hold onto it
    ping_count = 0;
    for np = 1:npings
        ping_num = burst(nb).ping(np); % ping number


        % read the corresponding binary file
        bin_file = sprintf('rec_7038_%06d.bin',ping_num);
        if ~exist(fullfile(swp_dir,bin_file),'file')
            fprintf('Warning: binary file missing for ping:%d\n',ping_num)
            continue
        end
        ip = find([png.num]==ping_num);
        % define sonar orientation angles
        yaw = (pi/180) * png(ip).rot_yaw;
        pitch = (pi/180)*png(ip).sen_pitch;
        roll = (pi/180)*png(ip).sen_roll;

        if(Verbose > 1)
            fprintf(' %s: pitch %0.2f, roll %0.2f, yaw %0.2f\n', bin_file, ...
                pitch*180/pi, roll*180/pi, yaw*180/pi);
        end

        % read raw element quadrature data
        try
            [hdr, data] = covis_read(fullfile(swp_dir, bin_file));
        catch
            fprintf('Warning: error reading ping %d at pitch %f\n',ping_num,burst(nb).pitch);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue
        end

        if isempty(data) || size(data,2)~=256
            fprintf('Warning: error reading ping %d at pitch %f\n',ping_num, burst(nb).pitch);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue;
        end

        if isempty(data)
            fprintf('Warning: error reading ping %d at pitch %f\n',ping_num, burst(nb).pitch);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue;
        end

        ping_count = ping_count+1;
        if(ping_count == 1)
            monitor = data;
            bf_sig_out = nan(size(data,1),size(data,2),npings);
        end

        % Correct phase using first ping as reference
        try
            data = covis_phase_correct(png(n), monitor, data);
        catch
            fprintf('Warning: error in phase correction for ping %d at pitch %f\n',ping_num,burst(nb).pitch);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue
        end

        % Apply Filter to data
        try
            [data, filt, png(n)] = covis_filter(data, filt, png(n));
        catch
            fprintf('Warning: error in filtering ping %d at pitch %f\n',ping_num,burst(nb).pitch);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue
        end

        % define beamformer parameters
        bfm.fc = png(n).hdr.xmit_freq;
        bfm.c = c;
        bfm.fs = png(n).hdr.sample_rate;
        bfm.first_samp = hdr.first_samp + 1;
        bfm.last_samp = hdr.last_samp;
        bfm.start_angle = -54;
        bfm.end_angle = 54;

        % beamform
        try
            [bfm, bf_sig] = covis_beamform(bfm, data);
        catch
            fprintf('Warning: error in beamforming for ping %d at pitch %f\n',ping_num,burst(nb).pitch);
            bad_ping_count = bad_ping_count+1;
            bad_ping(bad_ping_count) = ping_num;
            continue
        end
        bf_sig_out(:,:,np) = bf_sig;
    end

    if all(isnan(bf_sig_out(:)))
        fprintf('Warning: no valid pings at pitch %f\n',burst(nb).pitch);
        continue
    end

    burst_count = burst_count+1;
    if burst_count == 1
        xb_out = nan(bfm.num_beams,nbursts);
        yb_out = nan(size(xb_out));
        zb_out = nan(size(xb_out));
    end

    bf_sig_m = nanmean(bf_sig_out,3); % averaging over pings
    snr_m = 10*log10(abs(bf_sig_m).^2./noise_floor^2); % signal to noise ratio


    % OSCFAR detection section
    %magnitude squared of the beamformed output
    mag2 = bf_sig_m .* conj(bf_sig_m);
    % determine the specified quantile at each range (time) step
    clutter = prctile(mag2, clutterp, 2);
    dBclutter = 10*log10(clutter) * ones(1,length(bfm.angle));
    dBmag = 10*log10(mag2);
    dBSCR = dBmag - dBclutter; % signal to clutter ratio (SCR)
    indexD = dBSCR > scrthreshold; % detected samples
    bf_sig_m(~indexD)=sqrt(min(clutter));

    % calibration
    try
        bf_sig_m_cal = covis_calibration(bf_sig_m, bfm, png(n), cal,T, S, pH, lat,depth);
    catch
        fprintf('Warning: error in calibration at pitch %f\n',burst(nb).pitch);
        continue
    end

    % target scattering coefficient
    v = abs(bf_sig_m_cal).^2;

    % mask out data with low snr
    v(snr_m<snr_thresh) = min(v(:));

    % transform sonar coords into world coords
    range = bfm.range;
    azim = bfm.angle;

    % smooth the seafloor return in range
    v_filt = medfilt1(v,20);
    % loop over beams
    for m = 1:(bfm.num_beams)
        % find max return in each beam
        [~,ii] = max(v_filt(:,m));
        rb = range(ii);
        [xb_out(m,nb), yb_out(m,nb), zb_out(m,nb)] = covis_coords_darrell(origin,rb,azim(m),yaw,roll,pitch,central_head,pos.declination);
    end
end   % End loop over bursts

grd_out = covis.grid;
[grd_out.v,grd_out.w] = l2grid(xb_out,yb_out,zb_out,grd_out.x,grd_out.y,grd_out.v,grd_out.w);
n = find(grd_out.w);
grd_out.v(n) = grd_out.v(n)./grd_out.w(n);
grd_out.v(grd_out.v==0) = nan;


% save local copies of covis structs
covis_vers = covis_version();
covis.release = covis_vers.version_number;
covis.grid = grd_out;
covis.sweep = swp;
covis.ping = png;
covis.sonar.position = pos;
covis.processing.beamformer = bfm;
covis.processing.calibrate = cal;
covis.processing.filter = filt;
covis.processing.mask.noise_floor = noise_floor;
covis.processing.mask.snr = snr_thresh;
covis.processing.oscfar.clutterp = clutterp;
covis.processing.oscfar.scrthreshold = scrthreshold;
covis.burst = burst;
covis.bad_ping = bad_ping;


%% plot bathymetry
if fig == 1
    xg = covis.grid.x;
    yg = covis.grid.y;
    zb = covis.grid.v;
    figure
    pcolorjw(xg,yg,zb);
    axis image;
    caxis([-1 4]);
    xlabel('Easting of COVIS ( m )');
    ylabel('Northing of COVIS ( m )');
    h = colorbar;
    colormap('jet');
    title(h,'m');
end
