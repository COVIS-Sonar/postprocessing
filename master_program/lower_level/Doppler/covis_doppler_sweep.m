

path(path,'../../Common');

global Verbose;

% INITIALIZATION

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
snr_thresh = 50;

% OSCFAR parameters
clutterp = 65; % start with the 65% quantile
scrthreshold = 20; % dB threshold for signal-to-clutter ratio

% determine nominal rotator yaw
sy = str2double(swp_name(end));
switch(sy)
    case 1
        yaw0 = 135;
    case 2
        yaw0 = 70;
    case 3
        yaw0 = 210;
end


% START OF MAIN PROGRAM

% set up sweep directory
swp_dir = fullfile(swp_path, swp_name);

% check that archive dir exists
if(~exist(swp_dir,'dir'))
    error('Sweep directory does not exist\n');
end

% parse sweep.json file in data archive
swp_file = 'sweep.json';
if(~exist(fullfile(swp_dir, swp_file),'file'))
    error('sweep.json file does not exist\n');
end
json_str = fileread(fullfile(swp_dir, swp_file));
swp = parse_json(json_str);

% set sweep path and name
swp.path = swp_path;
swp.name = swp_name;

% parsing the json input file for the user supplied parameters
if(isempty(json_file) || (json_file == 0))
    % default json input file
    json_file = fullfile('..\..\input', 'covis_doppler.json');
end
% check that josn input file exists
if ~exist(json_file,'file')
    error('JSON input file does not exist')
end
% parse json file
json_str = fileread(json_file);
covis = parse_json(json_str);
if(strcmpi(covis.type, 'doppler') == 0)
    fprintf('Incorrect covis input file type\n');
    return;
end
Verbose = covis.user.verbose;


% define a 2D rectangular data grid
for n=1:length(covis.grid)
    [covis.grid{n}] = covis_rectgrid_doppler(covis.grid{n}); % corr grid
    covis.grid{n}.name = swp.name;
end

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

% Set compass declination
if(~isfield(pos,'declination'))
    pos.declination = 18.0;
end

% calibration parameters
if(~isfield(cal,'mode'))
    cal.mode = 'VSS'; % 'VSS' or 'TS'
end

% set position of sonar
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
if length(json_file)~=nfiles
    sprintf('there are %d json files missing in the sweep', nfiles-length(json_file));
end

% Read index file
% The index file contains the sweep parameters:
% ping,seconds,microseconds,pitch,roll,yaw,kPAngle,kRAngle,kHeading
% in csv format
ind_file = 'index.csv';
if ~exist(fullfile(swp_dir,ind_file),'file')
    error('Index file does not exists')
end
csv = csvread(fullfile(swp_dir, ind_file), 1, 0);
if(size(csv,1) ~= nfiles)
    fprintf('index size and file number mismatch');
end
if(Verbose)
    fprintf('Parsing %s\n', fullfile(swp_dir, ind_file));
end

% save index data in ping structure
central_yaw = 135;
central_head = 270;
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
    png(n).tcm.kHeading = (pos.declination + png(n).rot.yaw-central_yaw+central_head)';
end
% find the number of bursts in the sweep and number of pings in each burst
% by checking when the elevation changes
[ burst ] = covis_parse_bursts( png );

% range of elev steps to process (in degrees)
elev_start = (covis.processing.bounds.pitch.start);
elev_stop = (covis.processing.bounds.pitch.stop);


% loop over bursts
covarsum = 0;
nbursts = length(burst);
nb_count = 0;
for nb = 1:nbursts
    
    % check elevation
    if((burst(nb).elev < elev_start) || (burst(nb).elev > elev_stop))
        continue;
    end
    nb_count = nb_count+1;
    if(Verbose)
        fprintf('Burst %d: Elevation %f\n', nb, burst(nb).elev);
    end;
    
    % loop over pings in a burst
    bf_sig_out = zeros(0);
    for np = 1:burst(nb).npings
        
        % parse filenames
        n = burst(nb).start_ping + (np-1);
        bin_file = file(n).name;
        [type,ping_num] = strread(bin_file,'rec_%d_%d.bin');
        
        % read ping meta data from json file
        json_str = fileread(fullfile(swp_dir, json_file(1).name));
        json = parse_json(json_str);
        if abs(json.hdr.sample_rate-34482)>1000
            json_str = fileread(fullfile(swp_dir, json_file(2).name));
            json = parse_json(json_str);
        end
        png(n).hdr = json.hdr;
        
        % define sonar attitude (use TCM6 angles)
        elev = (pi/180) * png(n).tcm.kPAngle;
        roll = (pi/180) * png(n).tcm.kRAngle;
        yaw = (pi/180) * png(n).tcm.kHeading;
        
        if(Verbose > 1)
            fprintf(' %s: elev %f, roll %f, yaw %f\n', bin_file, ...
                elev*180/pi, roll*180/pi, yaw*180/pi);
        end;
        
        % read raw element quadrature data
        try
            [hdr, data] = covis_read(fullfile(swp_dir, bin_file));
        catch err
            disp(['bad ping: ',num2str(ping_num),' at elev angle: ',num2str(burst(nb).elev)])
            bf_sig_out(:,:,np) = nan;
            continue
        end
        
        if(np == 1)
            monitor = data;
        end
        
        if(strfind(dsp.ping_combination.mode,'diff'))
            % Correct phase using first ping as reference
            try
                data = covis_phase_correct(png(n), monitor, data);
            catch err
                disp(['bad ping: ',num2str(ping_num),' at elev angle: ',num2str(burst(nb).elev)])
                bf_sig_out(:,:,np) = nan;
                continue
            end
        end
        
        % Calculate the covariance function of the monitor
        % signal. The covariance function will be used to correct the
        % offset in the radial velocity estimates caused by the timing
        % mismatch between the clocks in COVIS for data generation and
        % digitization.
        try
            [covar] = covis_offset_covar(data,png(n));
        catch err
            disp(['bad ping: ',num2str(ping_num),' at elev angle: ',num2str(burst(nb).elev)])
            bf_sig_out(:,:,np) = nan;
        end
        covarsum = covarsum+covar;
        
        
        
        % Apply Filter to data
        [data, filt, png(n)] = covis_filter_doppler(data, filt, png(n));
        
        % define beamformer parameters
        bfm.fc = png(n).hdr.xmit_freq;
        bfm.c = c;
        bfm.fs = png(n).hdr.sample_rate;
        bfm.first_samp = hdr.first_samp + 1;
        bfm.last_samp = hdr.last_samp;
        yaw1 = pos.declination + yaw0-central_yaw+central_head;
        dang = yaw*180/pi-yaw1;
        bfm.start_angle = -54-dang;
        bfm.end_angle = 54-dang;
        
        % beamform
        try
            [bfm, bf_sig] = covis_beamform(bfm, data);
        catch err
            disp(['bad ping: ',num2str(ping_num),' at elev angle: ',num2str(burst(nb).elev)])
            bf_sig_out(:,:,np) = nan;
            continue
        end
        
        % calibration
        bf_sig_cal = covis_calibration(bf_sig, bfm, png(ip), cal,T, S, pH, lat,depth);
        
        bf_sig_out(:,:,np) = bf_sig_cal;
    end % End of loop over pings in burst
    
    if all(isnan(bf_sig_out))
        continue
    end
    
    % Compute Doppler shift using all pings within the burst
    range = bfm.range;
    [vr_cov,vr_vel,rc,I_av,vr_std,I_std,covar] = covis_incoher_dop_xgy(png(n).hdr, dsp, range, bf_sig_out,'diff');
    
    % Caculate snr and mask out data with low snr
    if nb_count==1
        bf_sig_noise = noise_floor*ones(size(bf_sig));
        bf_sig_noise = covis_calibration(bf_sig_noise,bfm,png(n),cal,T,S,pH,lat,depth);
        bf_sig_noise = repmat(bf_sig_noise,1,1,size(bf_sig_out,3));
        [~,~,~,I_noise,~,~,~] = covis_incoher_dop_xgy(png(n).hdr, dsp, range, bf_sig_noise,'ave');
    end
    snr = 10*log10(I_av./I_noise);
    I_av(snr<snr_thresh) = min(I_av(:));
    vr_cov(snr<snr_thresh) = nan;
    vr_vel(snr<snr_thresh) = nan;
    vr_std(snr<snr_thresh) = nan;
    I_std(snr<snr_thresh) = min(I_av(:));
    covar(snr<snr_thresh) = nan;
    
    
    % define sonar attitude
    elev = (pi/180) * png(n).tcm.kPAngle;
    roll = (pi/180) * png(n).tcm.kRAngle;
    yaw = (pi/180) * png(n).tcm.kHeading;
    
    % transform sonar coords into world coords
    azim = bfm.angle;
    [xv, yv, zv] = covis_coords(origin, rc, azim, yaw, roll, elev);
    
    % sin of elevation angle
    sin_elev = zv./sqrt(xv.^2+yv.^2+zv.^2);
    % Vertical velocity under assumption that flow is vertical
    vz_cov = vr_cov./sin_elev;
    
    % grid the data
    grd_in.x = xv;
    grd_in.y = yv;
    grd_in.z = zv;
    for n=1:length(covis.grid)
        grd_out = covis.grid{n};
        switch lower(grd_out.type)
            % grid the vertical velocity data
            case {'doppler velocity'}
                grd_in.vr_cov = vr_cov;
                grd_in.vr_vel = vr_vel;
                grd_in.vz_cov = vz_cov;
                grd_in.std = vr_std;
                grd_in.covar = covar;
                grd_in.offset_cov = covarsum;
                grd_out = l3grid_doppler(grd_in,grd_out);
                % grid the intensity data
            case {'intensity'}
                grd_in.I = I_av;
                grd_in.std = I_std;
                grd_out = l3grid_doppler(grd_in,grd_out);
            otherwise
                error('No Doppler grid is found')
        end
        covis.grid{n} = grd_out;
    end
    
end   % End of burst loop

% normalize the velocity grid with the grid weights
% save the grid structures in the main covis struct
for n=1:length(covis.grid)
    grd = covis.grid{n};
    switch lower(grd.type)
        case {'doppler velocity'}
            m = find(grd.w);
            grd.vr_cov(m) = grd.vr_cov(m)./grd.w(m);
            grd.vr_vel(m)=grd.vr_vel(m)./grd.w(m);
            grd.vz_cov(m)=grd.vz_cov(m)./grd.w(m);
            grd.std(m) = grd.std(m)./grd.w(m);
            grd.covar(m)=grd.covar(m)./grd.w(m);
            covis.grid{n} = grd;
        case {'intensity'}
            m = find(grd.w);
            grd.I(m) = grd.I(m)./grd.w(m);
            grd.std(m)=grd.std(m)./grd.w(m);
            covis.grid{n} = grd;
        otherwise
            error('No Doppler grid is found')
    end
end

% save the local copies of covis structs
covis.sweep = swp;
covis.ping = png;
covis.calibrate = pos;
covis.processing.beamformer = bfm;
covis.processing.calibrate = cal;
covis.processing.filter = filt;
covis.burst = burst;

% save results
out_dir = 'F:\COVIS\Axial\COVIS_data\Doppler\Sep\11\grid_data';
out_name = [swp_name,'_50dB_snr.mat'];
save(fullfile(out_dir,out_name),'covis');
%end






