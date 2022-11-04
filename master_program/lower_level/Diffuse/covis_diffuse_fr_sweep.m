% This function processes the raw data recorded in COVIS's
% bathymetry mode to generate gridded data for creating bathymetric and
% diffuse-flow maps


function covis = covis_diffuse_fr_sweep(swp_path, swp_name, json_file, fig)
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
% swp_path = 'F:\COVIS\Axial\COVIS_data\raw\full_range_imaging';
% swp_name = 'COVIS-20200628T000002-fullimaging1';
% json_file = 0;
% fig = 1;

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

% full path to sweep
swp_dir = fullfile(swp_path, swp_name);

% check that archive dir exists
if(~exist(swp_dir,'dir'))
    error('Sweep directory does not exist\n');
end

% parse sweep.json file in data archive
swp_file = 'sweep.json';
if(~exist(fullfile(swp_dir, swp_file),'file'))
    error('sweep.json file does not exist');
end
json_str = fileread(fullfile(swp_dir, swp_file));
swp = parse_json(json_str);

% set sweep path and name
swp.path = swp_path;
swp.name = swp_name;

% parsing the json input file for the user supplied parameters
if isempty(json_file) || all(json_file == 0)
    % default json input file
    json_file = find_input_file('covis_diffuse_fr.json');
end
% check that json input file exists
if ~exist(json_file,'file')
    error('JSON input file does not exist');
end

% parse json file
json_str = fileread(json_file);
covis = parse_json(json_str);
if(strcmpi(covis.type, 'diffuse_fr') == 0)
    error('Incorrect covis input file type. Looking for: diffuse_fr, Current type: %s\n',covis.type);
end

% Define global variable Verbose and read its value from the json file
global Verbose;
Verbose = covis.user.verbose;


% define a 2D rectangular data grid
for n=1:length(covis.grid)
    [covis.grid{n}] = covis_rectgrid(covis.grid{n}); % corr grid
    covis.grid{n}.name = swp.name;
end

% set local copies of covis structs
pos = covis.sonar.position;
bfm = covis.processing.beamformer;
cal = covis.processing.calibrate;
filt = covis.processing.filter;
cor = covis.processing.correlation;
avg = covis.processing.averaging;
noise = covis.processing.noise;



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
    cal.mode = 'TS-Fan'; % 'VSS' or 'TS'
end

% filter parameters
if(~isfield(cal,'filt'))
    cal.filt = 1;         % 1 for low-pass filtering, 0 for no filtering
    cal.filt_bw = 2.0;     % Bandwdith is filt_bw/tau
end

% set position of sonar
% should use covis.position info, but for now ...
if(~isfield(pos,'altitude'))
    pos.altitude = 4.2;
end
alt = pos.altitude;
origin = [0,0,alt]; % sonar is always at (0,0,0) in world coords

% correlation range window size [sec]
if(~isfield(cor,'window_size'))
    cor.window_size = 0.001;
end
window_size = cor.window_size;

% correlation window overlap [number of samples]
if(~isfield(cor,'window_overlap'))
    cor.window_overlap = 0.5;
end
window_overlap = cor.window_overlap;

% correlation time lag between pings
if(~isfield(cor,'tlag'))
    cor.tlag = 2;
end
tlag = cor.tlag;

% averaging window size (sec)
if(~isfield(avg,'avg_win'))
    avg.avg_win = 4;
end
avg_win = avg.avg_win;

% noise floor (in raw unit)
if(~isfield(noise,'noise_floor'))
    noise.noise_floor = 0.1970;
end
noise_floor = noise.noise_floor;

% Signal-to-noise threshold (dB)
if(~isfield(noise,'snr_thresh'))
    noise.snr_thresh = 90;
end
snr_thresh = noise.snr_thresh;

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

% numbers of samples in correlation window and window overlap
fsamp = json.hdr.sample_rate;
cwsize = round(fsamp*window_size);
cwovlap = round(window_overlap*cwsize);

% read index file
% The index file contains the sweep parameters:
% ping,seconds,microseconds,pitch,roll,yaw,kPAngle,kRAngle,kHeading
% in csv format
ind_file = 'index.csv';
if ~exist(fullfile(swp_dir,ind_file),'file')
    error('csv file does not exist in %', swp_dir)
else
    if(Verbose)
        fprintf('Parsing %s\n', fullfile(swp_dir, ind_file));
    end
    csv = csvread(fullfile(swp_dir, ind_file), 1, 0);
    csv = csv(csv(:,1)~=0,:); % remove the first dummy ping that logs sonar neutral orientation
end


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
        fprintf('Not enough pings in burst\n');
        continue;
    end
    
    % loop over pings in a burst, read data and hold onto it
    ping_count = 0;
    ping_num_burst = nan(1,npings);
    ping_sec_burst = nan(1,npings);
    for np = 1:npings
        ping_num = burst(nb).ping(np); % ping number
        ping_num_burst(np) = ping_num;
        ping_sec_burst(np) = png(ping_num).sec; % ping time (sec);
        % read the corresponding binary file
        bin_file = sprintf('rec_7038_%06d.bin',ping_num);
        if ~exist(fullfile(swp_dir,bin_file),'file')
            fprintf('binary file missing for ping:%d\n',ping_num)
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
            [hdr, data1] = covis_read(fullfile(swp_dir, bin_file));
        catch
            disp(['bad ping: ',num2str(ping_num),' at pitch: ',num2str(burst(nb).pitch)])
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue
        end
        
        if isempty(data1) || size(data1,2)~=256
            fprintf('Warning: error reading ping %d\n',ping_num);
            bad_ping_count = bad_ping_count + 1;
            bad_ping(bad_ping_count) = ping_num;
            continue;
        end
        
        ping_count = ping_count+1;
        if(ping_count == 1)
            data = nan(size(data1,1),size(data1,2),npings);
        end
        data(:,:,np) = data1;
    end  % End loop on pings
    
    if all(isnan(data(:)))
        fprintf('no valid pings at pitch %f\n',burst(nb).pitch);
        continue
    end
    
    % Loop again on pings to calculate backscatter intensity and
    % scintillation parameters
    I_t1 = nan(size(data));
    Isq_t = nan(size(data));
    bf_sig_t = nan(size(data));
    data0 = data(:,:,1);
    for np = 1:npings
        ping_num = ping_num_burst(np);
        if all(isnan(data0(:)))
            data0 = data(:,:,np+1);
            continue
        end
        data1 = data(:,:,np);
        % Correct phase
        try
            data1 = covis_phase_correct(png(ping_num), data0, data1);
        catch
            fprintf('Warning: error in phase correction for ping: %d\n',png(ping_num).num)
            continue
        end
        % define beamformer parameters
        bfm.fc = png(ping_num).hdr.xmit_freq;
        bfm.c = c;
        bfm.fs = png(ping_num).hdr.sample_rate;
        bfm.first_samp = hdr.first_samp + 1;
        bfm.last_samp = hdr.last_samp + 1;
        bfm.start_angle = -64;
        bfm.end_angle = 64;
        
        % Apply Filter to data
        try
            [data1, filt, png(ping_num)] = covis_filter(data1, filt, png(ping_num));
        catch
            fprintf('Warning: error in filtering ping: %d\n',png(ping_num).num)
            continue;
        end
        
        % beamform the quadrature data
        try
            [bfm, bf_sig1] = covis_beamform(bfm, data1);
        catch
            fprintf('Warning: error in beamforming for ping: %d\n',png(ping_num).num)
            continue;
        end
        
        % apply calibration to beamformed data
        try
            bf_sig1 = covis_calibration(bf_sig1, bfm, png(ping_num), cal,T,S,pH,lat,depth);
        catch
            fprintf('Warning: error in calibrating ping: %d\n',png(ping_num).num)
            continue;
        end
        
        bf_sig_t(:,:,np) = bf_sig1;
        
        % calculate scintillation index and log-amplitude flluctuations
        I1 = abs(bf_sig1).^2;
        Isq1 = abs(bf_sig1).^4;
        I_t1(:,:,np) = I1;
        Isq_t(:,:,np) = Isq1;
    end
    % average over the selected time window
    time_win1 = ping_sec_burst(1);
    time_win2 = time_win1+avg_win;
    ii_ave = find(ping_sec_burst>=time_win1&ping_sec_burst<=time_win2);
    nave = length(ii_ave); % number of samples used for averaging
    I_av1 = nanmean(I_t1(:,:,ii_ave),3);
    X = log(I_t1(:,:,ii_ave)./repmat(I_av1,1,1,length(ii_ave)));
    J = abs(nanmean(bf_sig_t(:,:,ii_ave)./abs(bf_sig_t(:,:,ii_ave)),3)).^2;
    X_var = nanvar(X(:,:,ii_ave),[],3);
    SI = nave*nansum(Isq_t(:,:,ii_ave),3)./nansum(I_t1(:,:,ii_ave),3).^2-1;
    Kp = nanmean(abs(bf_sig_t(:,:,ii_ave)).^2,3)./nanmean(abs(bf_sig_t(:,:,ii_ave)),3).^2-1;
    sig_phi2 = log(1./J);
    sig_phi2(J==0) = nan;
    
    % loop over pings again to form ping-ping decorrelation
    bf_sig_t_sub = bf_sig_t(:,:,ii_ave);
    ping_sec_sub = ping_sec_burst(ii_ave);
    ping_num_sub = ping_num_burst(ii_ave);
    ping_rate = png(ping_num_sub(1)).hdr.max_ping_rate;
    for np = 1:size(bf_sig_t_sub,3)
        ping_num = ping_num_sub(np);
        
        % get range and azimuthal angles
        range = bfm.range;
        azim = bfm.angle;
        
        % calculate the target strength corresponding to the noise floor
        if np==1
            bf_sig1 = squeeze(bf_sig_t(:,:,np));
            bf_sig_noise = noise_floor*ones(size(bf_sig1));
            bf_sig_noise = covis_calibration(bf_sig_noise,bfm,png(ping_num),cal,T,S,pH,lat,depth);
            [~, E1_noise,E2_noise, ~] = covis_covar_hamming(bf_sig_noise, bf_sig_noise, range, cwsize, cwovlap);
            I_noise1 = abs(bf_sig_noise).^2;
            I_noise2 = sqrt(E1_noise.*E2_noise);
            I_t2 = nan(size(I_noise2,1),size(I_noise2,2),size(data,3));
            cov_t = nan(size(I_t2));
            E1_t = nan(size(I_t2));
            E2_t = nan(size(I_t2));
        end
        
        t1 = png(ping_num).sec;
        t2 = t1+tlag;
        if ~isempty(find(abs(ping_sec_sub-t2)<1/ping_rate,1))
            [~,np2] = min(abs(ping_sec_sub-t2));
        else
            continue
        end
        
        % pair of pings used
        bf_sig1 = bf_sig_t_sub(:,:,np);
        bf_sig2 = bf_sig_t_sub(:,:,np2);
        
        % Correlate pings
        %  rc is the range of the center of the corr bin
        try
            [cov,E1,E2,rc] = covis_covar_hamming(bf_sig1, bf_sig2, range, cwsize, cwovlap);
        catch
            fprintf('Warning: error in calculating decorrelation between pings %d and %d\n',png(ping_num_sub(np)).num,png(ping_num_sub(np2)).num)
            continue
        end
        I2 = sqrt(E1.*E2);
        
        % save the quanities of interest
        cov_t(:,:,np) = cov;
        E1_t(:,:,np) = E1;
        E2_t(:,:,np) = E2;
        I_t2(:,:,np) = I2;
    end
    
    % average oveer ping pairs
    cor_av = abs(nansum(cov_t,3))./sqrt(nansum(E1_t,3).*nansum(E2_t,3));
    I_av2 = nanmean(I_t2,3);
    decor_I_av = (1-cor_av).*I_av2;
    
    % mask out data points with low SNR
    snr1 = 10*log10(I_av1./I_noise1);
    snr2 = 10*log10(I_av2./I_noise2);
    mask1 = nan(size(snr1));
    mask2 = nan(size(snr2));
    mask1(snr1>snr_thresh) = 1;
    mask2(snr2>snr_thresh) = 1;
    cor_av = cor_av.*mask2;
    decor_I_av = decor_I_av.*mask2;
    X_var = X_var.*mask1;
    Kp = Kp.*mask1;
    SI = SI.*mask1;
    sig_phi2 = sig_phi2.*mask1;
    
    
    burst_count = burst_count+1;
    if burst_count == 1
        xb_out = nan(bfm.num_beams,nbursts);
        yb_out = nan(size(xb_out));
        zb_out = nan(size(xb_out));
        xb2_out = nan(size(xb_out));
        yb2_out = nan(size(xb_out));
        zb2_out = nan(size(xb_out));
        I1_out = nan(size(xb_out));
        X_var_out = nan(size(xb_out));
        SI_out = nan(size(xb_out));
        Kp_out = nan(size(xb_out));
        sig_phi2_out = nan(size(xb_out));
        I2_out = nan(size(xb_out));
        cor_out = nan(size(xb_out));
        decor_I_out = nan(size(xb_out));
    end
    
    % % transform sonar coords into world coords
    
    % loop over beams
    for m = 1:(bfm.num_beams)
        % find max return in each beam
        [~,ii] = max(I_av1(:,m));
        rb = range(ii);
        I1_out(m,nb) = I_av1(ii,m);
        X_var_out(m,nb) = X_var(ii,m);
        SI_out(m,nb) = SI(ii,m);
        sig_phi2_out(m,nb) = sig_phi2(ii,m);
        Kp_out(m,nb) = Kp(ii,m);
        [xb_out(m,nb), yb_out(m,nb), zb_out(m,nb)] = covis_coords_darrell(origin,rb,azim(m),yaw,roll,pitch,central_head,pos.declination);
        [~,ii] = max(I_av2(:,m));
        rb2 = rc(ii);
        I2_out(m,nb) = I_av2(ii,m);
        cor_out(m,nb) = cor_av(ii,m);
        decor_I_out(m,nb) = decor_I_av(ii,m);
        [xb2_out(m,nb), yb2_out(m,nb), zb2_out(m,nb)] = covis_coords_darrell(origin,rb2,azim(m),yaw,roll,pitch,central_head,pos.declination);
    end
end   % End loop over bursts


% interpolate onto a 2D grid
for k = 1:length(covis.grid)
    grd = covis.grid{k};
    switch grd.type
        case 'bathy_filt'
            [grd.v, grd.w] = l2grid(xb2_out, yb2_out, zb2_out, grd.x, grd.y, grd.v, grd.w);
        case 'intensity_filt'
            [grd.v, grd.w] = l2grid(xb2_out, yb_out, I2_out, grd.x, grd.y, grd.v, grd.w);
        case 'decorrelation intensity'
            [grd.v, grd.w] = l2grid(xb2_out, yb2_out, decor_I_out, grd.x, grd.y, grd.v, grd.w);
        case 'decorrelation'
            [grd.v, grd.w] = l2grid(xb2_out, yb2_out, 1-cor_out, grd.x, grd.y, grd.v, grd.w);
        case 'bathy'
            [grd.v, grd.w] = l2grid(xb_out, yb_out, zb_out, grd.x, grd.y, grd.v, grd.w);
        case 'intensity'
            [grd.v, grd.w] = l2grid(xb_out, yb_out, I1_out, grd.x, grd.y, grd.v, grd.w);
        case 'Chi_var'
            [grd.v, grd.w] = l2grid(xb_out, yb_out, X_var_out, grd.x, grd.y, grd.v, grd.w);
        case 'SI'
            [grd.v, grd.w] = l2grid(xb_out, yb_out, SI_out, grd.x, grd.y, grd.v, grd.w);
        case 'Sig_phi2'
            [grd.v, grd.w] = l2grid(xb_out, yb_out, sig_phi2_out, grd.x, grd.y, grd.v, grd.w);
        case 'Kp'
            [grd.v, grd.w] = l2grid(xb_out, yb_out, Kp_out, grd.x, grd.y, grd.v, grd.w);
    end
    covis.grid{k} =grd;
end

% normalize the grid with the grid weights
for k = 1:length(covis.grid)
    grd = covis.grid{k};
    ii = find(grd.w);
    grd.v(ii) = grd.v(ii)./grd.w(ii);
    covis.grid{k} = grd;
end

% save metadata into the covis structure
covis_vers = covis_version();
covis.release = covis_vers.version_number;
covis.sweep = swp;
covis.ping = png;
covis.sonar.position = pos;
covis.processing.beamformer = bfm;
covis.processing.calibrate = cal;
covis.processing.filter = filt;
covis.processing.correlation.tlag = tlag;
covis.processing.averaging.win = avg_win;
covis.processing.snr.noise_floor = noise_floor;
covis.processing.snr.threshold = snr_thresh;
covis.bad_ping = bad_ping;
fclose('all');

%% plot results
xg = covis.grid{1}.x;
yg = covis.grid{1}.y;
if fig==1
    for k = 1:length(covis.grid)
        grd = covis.grid{k};
        switch grd.type
            case 'bathy'
                figure
                grd.v(grd.v==0) = nan;
                pcolorjw(xg,yg,grd.v);
                shading flat
                axis image;
                caxis([-1 4]);
                colormap('jet');
                h = colorbar;
                title(h,'m');
                xlabel('Easting of COVIS ( m )');
                ylabel('Northing of COVIS ( m )');
                hold on;
                plot(0,0,'.m','markersize',30);
                hold off;
                title('bathy');
            case 'intensity'
                figure
                grd.v(grd.v==0) = nan;
                pcolorjw(xg,yg,10*log10(grd.v));
                shading flat
                axis image;
                caxis([-50 -15]);
                colormap('jet');
                h = colorbar;
                title(h,'dB');
                xlabel('Easting of COVIS ( m )');
                ylabel('Northing of COVIS ( m )');
                hold on;
                plot(0,0,'.m','markersize',30);
                hold off;
                title('TS');
            case 'Kp'
                grd.v(grd.v==0) = nan;
                figure
                pcolorjw(xg,yg,grd.v);
                shading flat
                axis image;
                caxis([0 0.5]);
                colormap('jet');
                h = colorbar;
                xlabel('Easting of COVIS ( m )');
                ylabel('Northing of COVIS ( m )');
                hold on;
                plot(0,0,'.m','markersize',30);
                hold off;
                title('Kp');
                
             case 'bathy_filt'
                figure
                grd.v(grd.v==0) = nan;
                pcolorjw(xg,yg,grd.v);
                shading flat
                axis image;
                caxis([-1 4]);
                colormap('jet');
                h = colorbar;
                title(h,'m');
                xlabel('Easting of COVIS ( m )');
                ylabel('Northing of COVIS ( m )');
                hold on;
                plot(0,0,'.m','markersize',30);
                hold off;
                title('bathy filt');
                
            case 'intensity_filt'
                figure
                grd.v(grd.v==0) = nan;
                pcolorjw(xg,yg,10*log10(grd.v));
                shading flat
                axis image;
                caxis([-50 -15]);
                colormap('jet');
                h = colorbar;
                title(h,'dB');
                xlabel('Easting of COVIS ( m )');
                ylabel('Northing of COVIS ( m )');
                hold on;
                plot(0,0,'.m','markersize',30);
                hold off;
                title('TS_filt');
                
            case 'decorrelation'
                grd.v(grd.v==0) = nan;
                figure
                pcolorjw(xg,yg,grd.v);
                shading flat
                axis image;
                caxis([0 0.5]);
                colormap('jet');
                h = colorbar;
                xlabel('Easting of COVIS ( m )');
                ylabel('Northing of COVIS ( m )');
                hold on;
                plot(0,0,'.m','markersize',30);
                hold off;
                title('decorrelation');
                
        end
    end
end
end
    