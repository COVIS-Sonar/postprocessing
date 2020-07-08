% This function process the raw data recorded in COVIS's Imaging mode to
% generate gridded data that is saved in the 'covis' structure array

% version 1.0 by guangyux@uw.edu (Oct 19, 2019)
%  --based on the original code written by Chris Jones in 2010

%function covis = covis_imaging_sweep(swp_path, swp_name, json_file, fig)
% Input:
% swp_path: raw data directory
% swp_name: name of raw data sweep
% json_file: full path and name of the json file that includes the metadata needed for
% processing. if 0, use the default json file.
% fig: set fig to 1 to plot gridded data. Plotting is muted otherwise.
% Output:
% covis: this is the Matlab structure array that includes the gridded data
% and metadata

% Example
swp_path = 'F:\COVIS\Axial\COVIS_data\raw\Doppler';
swp_name = 'COVIS-20191113T020002-doppler1';
json_file = 0;



%% Initialization

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


% snr parameters
noise_floor = 0.1970; % rms noise floor
snr_thresh = 50; % threshold (dB)


% OSCFAR parameters
clutterp = 65; % start with the 65% quantile as the clutter threshold
scrthreshold = 20; % dB threshold for signal-to-clutter ratio


% full path to sweep
swp_dir = fullfile(swp_path, swp_name);

% check that archive dir exists
if(~exist(swp_dir,'dir'))
    error('Sweep directory does not exist\n');
end

% parse sweep.json file in data archive
swp_file = 'sweep.json';
if(~exist(fullfile(swp_dir, swp_file),'file'))
    disp('sweep.json file does not exist');
    return;
end
json_str = fileread(fullfile(swp_dir, swp_file));
swp = parse_json(json_str);

% save sweep path and name in swp structure
swp.path = swp_path;
swp.name = swp_name;

% parsing the json input file for the user supplied parameters
if isempty(json_file) || all(json_file == 0)
    % default json input file
    json_file = find_input_file('covis_doppler.json');
end
% check that json input file exists
if ~exist(json_file,'file')
    error('JSON input file does not exist');
end
% parse json file
json_str = fileread(json_file);
covis = parse_json(json_str);
if(strcmpi(covis.type, 'doppler') == 0)
    fprintf('Incorrect covis input file type\n');
    return;
end

% Define global variable Verbose and read its value from the json file
global Verbose;
Verbose = covis.user.verbose;

% define a 3D rectangular data grid
for n=1:length(covis.grid)
    [covis.grid{n}] = covis_rectgrid_doppler(covis.grid{n}); % corr grid
    covis.grid{n}.name = swp.name;
end
grd_out = covis.grid;

% set local copies of covis structs
grd = covis.grid;
usr = covis.user;
pos = covis.sonar.position;
dsp = covis.processing;
bfm = covis.processing.beamformer;
cal = covis.processing.calibrate;
filt = covis.processing.filter;



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
    fprintf('there are %d json files missing in the sweep', nfiles-length(json_file));
elseif length(json_file)>nfiles
    fprintf('there are %d bin files missing in the sweep', length(json_file)-nfiles);
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
csv = csv(csv(:,1)~=0,:); % remove the first dummy ping that logs sonar neutral orientations
% save index data in png structure
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
        fprintf('Burst %d: pitch %f\n', nb, burst(nb).pitch);
    end
    
    npings = burst(nb).npings;
    
    % check that there's enough pings in burst
    if((npings < 2))
        fprintf('Not enough pings in burst\n');
        continue;
    end
    
    % loop over pings in a burst
    ping_count = 0;
    for np = 1:npings
        ping_num = burst(nb).ping(np); % ping number
        % read the corresponding binary file
        bin_file = sprintf('rec_7038_%06d.bin',ping_num);
        if ~exist(fullfile(swp_dir,bin_file),'file')
            fprintf('binary file missing for ping:%d\n',ping_num)
            continue
        end
        ip = find([png.num]==ping_num);
        % define sonar orientation based on TCM readings
        pitch = (pi/180) * png(ip).sen_pitch;
        roll = (pi/180) * png(ip).sen_roll;
        yaw = (pi/180) * png(ip).rot_yaw;
        
        if(Verbose > 1)
            fprintf('Reading %s\n', fullfile(swp_dir, bin_file));
        end
        if(Verbose > 2)
            fprintf('Pitch %f, roll %f, yaw %f\n',pitch*180/pi, roll*180/pi, yaw*180/pi);
        end
        
        % read raw element quadrature data
        try
            [hdr, data] = covis_read(fullfile(swp_dir, bin_file));
        catch
            fprintf('error reading ping %d at pitch %f\n',ping_num,burst(nb).pitch);
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
        
        ping_count = ping_count+1;
        if(ping_count == 1)
            monitor = data;
            bf_sig_out = nan(size(data,1),size(data,2),npings);
        end
        
        % Correct phase using first ping as reference
        try
            data = covis_phase_correct(png(ip), monitor, data);
        catch
            fprintf('error in phase correction for ping %d at pitch %f\n',ping_num,burst(nb).pitch);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue
        end
        
        % Apply Filter to data
        try
            [data, filt, png(ip)] = covis_filter(data, filt, png(ip));
        catch
            fprintf('error in filtering ping %d at pitch %f\n',ping_num,burst(nb).pitch);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue
        end
        
        % define beamformer parameters
        bfm.fc = png(ip).hdr.xmit_freq;
        bfm.c = c;
        bfm.fs = png(ip).hdr.sample_rate;
        bfm.first_samp = hdr.first_samp + 1;
        bfm.last_samp = hdr.last_samp;
        bfm.start_angle = -54;
        bfm.end_angle = 54;
        
        % conduct beamforming
        try
            [bfm, bf_sig] = covis_beamform(bfm, data);
        catch
            fprintf('error in beamforming for ping %d at pitch %f\n',ping_num,burst(nb).pitch);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue
        end
        % calibration
        try
            bf_sig = covis_calibration(bf_sig, bfm, png(n), cal,T, S, pH, lat,depth);
        catch
            fprintf('Warning: error in calibration at pitch %f\n',burst(nb).pitch);
            continue
        end
            
        bf_sig_out(:,:,np) = bf_sig;
    end
      
    if all(isnan(bf_sig_out(:)))
        fprintf('no valid pings at pitch %f\n',burst(nb).pitch);
        continue
    end
    
    burst_count = burst_count+1;
    if burst_count == 1
        
        xv_out1 = nan(size(bf_sig,1),size(bf_sig,2),nbursts);
        yv_out1 = nan(size(xv_out1));
        zv_out1 = nan(size(xv_out1));
        Id_out1 = nan(size(xv_out1));
        Id_filt_out1 = nan(size(xv_out1));
        
        
        cor = dsp.correlation;
        % correlation range window size [number of samples]
        if(~isfield(cor,'window_size'))
            cor.window_size = 0.001;
        end
        window_size = cor.window_size;
        
        % correlation window overlap [number of samples]
        if(~isfield(cor,'window_overlap'))
            cor.window_overlap = 0.5;
        end
        overlap = cor.window_overlap;
        dsp.cor = cor;
        
        fsamp = png(ip).hdr.sample_rate;
        nwindow = round(window_size*fsamp);
        noverlap = round(overlap*nwindow);
        nbins=floor(size(bf_sig_out,1)/(nwindow-noverlap));
        xv_out2 = nan(nbins,size(bf_sig_out,2),nbursts);
        yv_out2 = nan(size(xv_out2));
        zv_out2 = nan(size(xv_out2));
        vr_cov_out = nan(size(xv_out2));
        vr_vel_out = nan(size(xv_out2));
        Id_out2 = nan(size(xv_out2));
        Id_filt_out2 = nan(size(xv_out2));
        vr_std_out = nan(size(xv_out2));
        covar_out = nan(size(xv_out2));
    end
    

    
    % ping average
    average = nanmean(bf_sig_out,3);
    bf_sig_d = sqrt(nanmean(abs(bf_sig_out-repmat(average,1,1,size(bf_sig_out,3))).^2,3)); % remove the average to enhance plume signals
    
    
    % signal-to-noise ratio
    noise_sig = noise_floor*ones(size(bf_sig_d));
    noise_sig = covis_calibration(noise_sig, bfm, png(n), cal,T, S, pH, lat,depth);
    snr_d = 20*log10(abs(bf_sig_d)./noise_sig);
     
    % OSCFAR detection section
    mag2 = bf_sig_d .* conj(bf_sig_d); % magnitude squared of the beamformed output
    % determine the specified quantile at each range (time) step
    clutter = prctile(mag2, clutterp, 2);
    dBclutter = 10*log10(clutter) * ones(1,length(bfm.angle));
    dBmag = 10*log10(mag2);
    dBSCR = dBmag - dBclutter; % signal to clutter ratio (SCR)
    indexD_d = dBSCR > scrthreshold; % detected samples
    

    
    Id = abs(bf_sig_d).^2;
    Id_filt = Id;
    
    % mask out data with low snr
    Id_filt(snr_d<snr_thresh) = 10^-9;
    Id_filt(~indexD_d) = 10^-9;
    Id(snr_d<snr_thresh) = 10^-9;

    
    % Compute Doppler shift using all pings within the burst
    range = bfm.range;
    [vr_cov,vr_vel,rc,Id2,vr_std,covar] = covis_incoher_dop_xgy(png(ip).hdr, dsp, range, bf_sig_out);
    Id_filt2 = Id2;
        
    % calculate snr and mask out backscatter with low snr
    [~,~,~,I_noise] = covis_incoher_dop_xgy(png(ip).hdr, dsp, range, noise_sig);
    snr_d2 = 10*log10(Id2./I_noise);

    % OSCFAR detection section
    mag2 = Id2; % magnitude squared of the beamformed output
    % determine the specified quantile at each range (time) step
    clutter = prctile(mag2, clutterp, 2);
    dBclutter = 10*log10(clutter) * ones(1,length(bfm.angle));
    dBmag = 10*log10(mag2);
    dBSCR = dBmag - dBclutter; % signal to clutter ratio (SCR)
    indexD_d2 = dBSCR > scrthreshold; % detected samples
    
    % mask out data with low snr
    Id_filt2(snr_d2<snr_thresh) = 10^-9;
    Id_filt2(~indexD_d2) = 10^-9;
    Id2(snr_d2<snr_thresh) = 10^-9;
    vr_cov(snr_d2<snr_thresh) = nan;
    vr_vel(snr_d2<snr_thresh) = nan;
    vr_std(snr_d2<snr_thresh) = nan;
    covar(snr_d2<snr_thresh) = nan;
    
    % transform sonar coords into world coords
    azim = bfm.angle;
    [xv, yv, zv] = covis_coords_darrell(origin, range, azim, yaw, roll, pitch, central_head, pos.declination);
    xv_out1(:,:,nb) = xv;
    yv_out1(:,:,nb) = yv;
    zv_out1(:,:,nb) = zv;
    Id_out1(:,:,nb) = Id;
    Id_filt_out1(:,:,nb) = Id_filt;
    
    [xv2, yv2, zv2] = covis_coords_darrell(origin, rc, azim, yaw, roll, pitch, central_head, pos.declination);
    xv_out2(:,:,nb) = xv2;
    yv_out2(:,:,nb) = yv2;
    zv_out2(:,:,nb) = zv2;
    Id_out2(:,:,nb) = Id2;
    Id_filt_out2(:,:,nb) = Id_filt2;
    vr_cov_out(:,:,nb) = vr_cov;
    vr_vel_out(:,:,nb) = vr_vel;
    vr_std_out(:,:,nb) = vr_std;
    covar_out(:,:,nb) = covar;
end   % End loop over bursts


% grid the data

for n=1:length(covis.grid)
    grd_out = covis.grid{n};
    switch lower(grd_out.type)
        % grid the vertical velocity data
        case {'doppler velocity'}
            grd_in.x = xv_out2;
            grd_in.y = yv_out2;
            grd_in.z = zv_out2;
            grd_in.vr_cov = vr_cov_out;
            grd_in.vr_vel = vr_vel_out;
            grd_in.std = vr_std_out;
            grd_in.covar = covar_out;
            grd_out = l3grid_doppler(grd_in,grd_out);
            % grid the intensity data
        case {'intensity'}
            grd_in.x = xv_out1;
            grd_in.y = yv_out1;
            grd_in.z = zv_out1;
            grd_in.Id = Id_out1;
            grd_in.Id_filt = Id_filt_out1;
            grd_out = l3grid_doppler(grd_in,grd_out);          
        case {'intensity_win'}
            grd_in.x = xv_out2;
            grd_in.y = yv_out2;
            grd_in.z = zv_out2;
            grd_in.Id = Id_out2;
            grd_in.Id_filt = Id_filt_out2;
            grd_out = l3grid_doppler(grd_in,grd_out);
        otherwise
            error('No Doppler grid is found')
    end
    covis.grid{n} = grd_out;
end

% normalize the velocity grid with the grid weights
for n=1:length(covis.grid)
    grd = covis.grid{n};
    switch lower(grd.type)
        case {'doppler velocity'}
            m = find(grd.w);
            grd.vr_cov(m) = grd.vr_cov(m)./grd.w(m);
            grd.vr_vel(m)=grd.vr_vel(m)./grd.w(m);
            grd.std(m) = grd.std(m)./grd.w(m);
            grd.covar(m)=grd.covar(m)./grd.w(m);
            covis.grid{n} = grd;
        case {'intensity'}
            m = find(grd.w);
            grd.Id(m) = grd.Id(m)./grd.w(m);
            grd.Id_filt(m)=grd.Id_filt(m)./grd.w(m);
            covis.grid{n} = grd;  
        case {'intensity_win'}
            m = find(grd.w);
            grd.Id(m) = grd.Id(m)./grd.w(m);
            grd.Id_filt(m)=grd.Id_filt(m)./grd.w(m);
            covis.grid{n} = grd;
        otherwise
            error('No Doppler grid is found')
    end
end

% save local copies of covis structs
covis_vers = covis_version();
covis.release = covis_vers.version_number;
covis.sweep = swp;
covis.ping = png;
covis.sonar.position = pos;
covis.processing.beamformer = bfm;
covis.processing.calibrate = cal;
covis.processing.filter = filt;
covis.processing.snr.noise_floor = noise_floor;
covis.processing.snr.threshold = snr_thresh;
covis.processing.oscfar.clutterp = clutterp;
covis.processing.oscfar.scrthreshold = scrthreshold;
covis.burst = burst;
covis.bad_ping = bad_ping;


