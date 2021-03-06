% This function process the raw data recorded in COVIS's Diffuse-flow mode to
% generate gridded data that is saved in the 'covis' structure array

% version 1.0 by guangyux@uw.edu (Oct 19, 2019)
%  --based on the original code written by Chris Jones in 2010
% version 2.0 by guangyux@uw.edu (Feb 27, 2020)
%  --use fixed time window to calculate averages
%  --use fixed time lag to calculate correlation
%  --more adjustable parameters used in processing are added to the metadata in the output structure
%  --version number of the code is added to the metadata in the output structure

function covis = covis_diffuse_sweep(swp_path, swp_name, json_file, fig)
% Input
%  swp_path: sweep directory
%  swp_name: sweep name
%  json_file: json file directory (if 0, use default json file)
%  fig: set fig to 1 to plot gridded data. plotting is muted otherwise.
% Output
%  covis: Matlab structure array that contains the gridded data and metadata

% Example
% swp_path = 'C:\COVIS\Axial\COVIS_data\raw\raw_data_combine\2019\08\01';
% swp_name = 'COVIS-20190801T003002-diffuse1';
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

% mask parameters
noise_floor = 0.64; % rms noise floor (uncalibrated in machine units)
snr_thresh = 45; % snr threshold ( dB )


% correlation time lag (sec)
tlag = 2;

% averaging window length (sec)
avg_win = 4;



% bathymetry data directory
swp_date = datenum(swp_name(7:21),'yyyymmddTHHMMSS');
if swp_date<=datenum(2019,7,6)
    bathy_file = sprintf('covis_bathy_2018.mat');
elseif swp_date<=datenum(2019,11,23)
    bathy_file = sprintf('covis_bathy_2019a.mat');
else
    bathy_file = sprintf('covis_bathy_2019b.mat');
end
bathy = load(find_input_file(bathy_file));




% full path to sweep
swp_dir = fullfile(swp_path, swp_name);

% check that archive dir exists
if(~exist(swp_dir,'dir'))
    error('Sweep directory %s does not exist\n',swp_dir);
end


% parse sweep.json file in data archive
swp_file = 'sweep.json';
if(~exist(fullfile(swp_dir, swp_file),'file'))
    error('sweep.json file does not exist in %s\n', swp_dir);
end
json_str = fileread(fullfile(swp_dir, swp_file));
swp = parse_json(json_str);

% save sweep path and name in swp structure
swp.path = swp_path;
swp.name = swp_name;

% parsing the json input file for the user supplied parameters
if(isempty(json_file) || all(json_file == 0))
    % default json input file
    json_file = find_input_file('covis_diffuse.json');
end
% check that json input file exists
if ~exist(json_file,'file')
    error('JSON input file does not exist');
end
json_str = fileread(json_file);
covis = parse_json(json_str);
if(strcmpi(covis.type, 'diffuse') == 0)
    error('Incorrect covis input file type. Looking for: diffuse, Current type: %s\n',covis.type);
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
    cal.mode = 'TS-Wide'; % 'VSS', 'TS-Wide', 'TS-Fan'
end
if(~isfield(cal,'filt'))
    cal.filt = 1;         % 1 for low-pass filtering, 0 for no filtering
    cal.filt_bw = 2.0;     % Bandwdith is filt_bw/tau
end

% sonar height off bottom
if(~isfield(pos,'altitude'))
    pos.altitude = 4.2;
end

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


% directory list of *.bin file
file = dir(fullfile(swp_dir, '*.bin'));
nfiles = length(file);

% check if there are missing json or binary files
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
while abs(json.hdr.sample_rate-34482)>1 && ind<=length(json_file)
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

% sonar orientation during the sweep
pitch = (pi/180) * png(1).sen_pitch;
roll = (pi/180) * png(1).sen_roll;
yaw = (pi/180) * png(1).rot_yaw;


%% Main program
% Loop over ping files
bad_ping = zeros(0);
bad_ping_count = 0;
count = 0;
nping = size(csv,1);
for np = 1:nping
    ping_num = png(np).num;
    bin_file = sprintf('rec_7038_%06d.bin',ping_num);
    if ~exist(fullfile(swp_dir,bin_file),'file')
        fprintf('Warning: binary file missing for ping:%d\n',ping_num)
        continue
    end

    if (np == 1)
        if(Verbose > 1)
            png(np).hdr  % View essential parameters
        end
    end

    if(Verbose > 1)
        fprintf('Reading %s: pitch %f, roll %f, yaw %f\n', bin_file, ...
            pitch*180/pi, roll*180/pi, yaw*180/pi);
    end

    % read raw element quadrature data
    try
        [hdr, data1] = covis_read(fullfile(swp_dir, bin_file));
    catch
        fprintf('Warning: error reading ping %d\n',ping_num);
        bad_ping_count = bad_ping_count + 1;
        bad_ping(bad_ping_count) = ping_num;
        continue;
    end
    if isempty(data1) || size(data1,2)~=128
        fprintf('Warning: error reading ping %d\n',ping_num);
        bad_ping_count = bad_ping_count + 1;
        bad_ping(bad_ping_count) = ping_num;
        continue;
    end
    count = count+1;
    if count == 1
        data = nan(size(data1,1),size(data1,2),nping);
    end
    data(:,:,np) = data1;
end   % End loop on pings


% Loop again on pings to calculate backscatter intensity and
% scintillation parameters
I_t1 = nan(size(data));
Isq_t = nan(size(data));
bf_sig_t = nan(size(data));
data0 = data(:,:,1);
for np = 1:size(data,3)
    if all(isnan(data0(:)))
        data0 = data(:,:,np+1);
        continue
    end
    data1 = data(:,:,np);
    % Correct phase
    try
        data1 = covis_phase_correct(png(np), data0, data1);
    catch
        fprintf('Warning: error in phase correction for ping: %d\n',png(np).num)
        continue
    end
    % define beamformer parameters
    bfm.fc = png(np).hdr.xmit_freq;
    bfm.c = c;
    bfm.fs = png(np).hdr.sample_rate;
    bfm.first_samp = hdr.first_samp + 1;
    bfm.last_samp = hdr.last_samp + 1;
    bfm.start_angle = -64;
    bfm.end_angle = 64;

    % Apply Filter to data
    try
        [data1, filt, png(np)] = covis_filter(data1, filt, png(np));
    catch
        fprintf('Warning: error in filtering ping: %d\n',png(np).num)
        continue;
    end

    % beamform the quadrature data
    try
        [bfm, bf_sig1] = covis_beamform(bfm, data1);
    catch
        fprintf('Warning: error in beamforming for ping: %d\n',png(np).num)
        continue;
    end
    % apply calibration to beamformed data
    try
        bf_sig1 = covis_calibration(bf_sig1, bfm, png(np), cal,T,S,pH,lat,depth);
    catch
        fprintf('Warning: error in calibrating ping: %d\n',png(np).num)
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
ping_sec = [png(1:nping).sec];
time_win1 = ping_sec(1);
time_win2 = time_win1+avg_win;
ii_ave = find(ping_sec>=time_win1&ping_sec<=time_win2);
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
ping_rate = png(1).hdr.max_ping_rate;
bf_sig_t_sub = bf_sig_t(:,:,ii_ave);
ping_sec_sub = ping_sec(ii_ave);
for np = 1:size(bf_sig_t_sub,3)


    % get range and azimuthal angles
    range = bfm.range;
    azim = bfm.angle;

    % calculate the target strength corresponding to the noise floor
    if np==1
        bf_sig1 = squeeze(bf_sig_t(:,:,np));
        bf_sig_noise = noise_floor*ones(size(bf_sig1));
        bf_sig_noise = covis_calibration(bf_sig_noise,bfm,png(n),cal,T,S,pH,lat,depth);
        [~, E1_noise,E2_noise, ~] = covis_covar_hamming(bf_sig_noise, bf_sig_noise, range, cwsize, cwovlap);
        I_noise1 = abs(bf_sig_noise).^2;
        I_noise2 = sqrt(E1_noise.*E2_noise);
        I_t2 = nan(size(I_noise2,1),size(I_noise2,2),size(data,3));
        cov_t = nan(size(I_t2));
        E1_t = nan(size(I_t2));
        E2_t = nan(size(I_t2));
    end

    t1 = png(np).sec;
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
        fprintf('Warning: error in calculating decorrelation between pings %d and %d\n',png(np).num,png(np2).num)
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



%% Find seafloor origins of acoustic backscatter

% load bathymetry data
alt = covis.sonar.position.altitude; % height of COVIS
covis_bathy = bathy.covis;
x_bathy = covis_bathy.grid.x;
y_bathy = covis_bathy.grid.y;
z_bathy = covis_bathy.grid.v-alt;



% Determine the real-world coordinates of pesudo planes perpendicular along each
% beam and the intercept of the plane on the seafloor

ver_beam = (-40:0.5:40)*pi/180; % vertical beam width

% rotation matrices
R = [cos(roll) 0 -sin(roll) ; 0 1 0 ; sin(roll) 0 cos(roll)];
Y = [cos(yaw) -sin(yaw) 0 ; sin(yaw) cos(yaw) 0 ; 0 0 1];
P = [1 0 0 ; 0 cos(pitch) sin(pitch) ;  0 -sin(pitch) cos(pitch)];

% COVIS central pointing direction re true north
ang = pi/180*(central_head + pos.declination -360);
ROT = [cos(ang) sin(ang) 0 ; -sin(ang) cos(ang) 0 ; 0 0 1];
M = ROT*Y'*P'*R';
x_out1 = zeros(length(range),length(azim));
y_out1 = zeros(length(range),length(azim));
x_out2 = zeros(length(rc),length(azim));
y_out2 = zeros(length(rc),length(azim));

for j = 1:length(azim)

    % receiver coordinates of a pseudo plane perpendicular along the beam
    azim1 = azim(j);
    xr = range(:)*cos(ver_beam)*sin(azim1);
    yr = range(:)*cos(ver_beam)*cos(azim1);
    zr = range(:)*sin(ver_beam);
    rr = [xr(:)';yr(:)';zr(:)'];

    % real-world coordinates of a pseudo plane
    rw = M*rr;

    xw1 = reshape(rw(1,:),length(range),length(ver_beam));
    yw1 = reshape(rw(2,:),length(range),length(ver_beam));
    zw1 = reshape(rw(3,:),length(range),length(ver_beam));

    % Find the intercept of the pseudo plane on the seafloor
    zw2 = interp2(x_bathy,y_bathy,z_bathy,xw1,yw1,'linear',nan);
    dz = abs(zw1-zw2);
    ii = find(dz<0.1);
    xw_sub = xw1(ii);
    yw_sub = yw1(ii);
    zw2_sub = zw2(ii);
    L_sub = sqrt(xw_sub.^2+yw_sub.^2+zw2_sub.^2);
    for i = 1:length(range)
        range1 = range(i);
        ii = find(abs(L_sub-range1)<0.5*min(diff(rc)));
        if ~isempty(ii)
            x_out1(i,j) = mean(xw_sub(ii));
            y_out1(i,j) = mean(yw_sub(ii));
        else
            x_out1(i,j) = nan;
            y_out1(i,j) = nan;
        end
    end
    for i = 1:length(rc)
        rc1 = rc(i);
        ii = find(abs(L_sub-rc1)<0.5*min(diff(rc)));
        if ~isempty(ii)
            x_out2(i,j) = mean(xw_sub(ii));
            y_out2(i,j) = mean(yw_sub(ii));
        else
            x_out2(i,j) = nan;
            y_out2(i,j) = nan;
        end
    end
end

% interpolate onto a 2D grid
for k = 1:length(covis.grid)
    grd = covis.grid{k};
    switch grd.type
        case 'decorrelation intensity'
            [grd.v, grd.w] = l2grid(x_out2, y_out2, decor_I_av, grd.x, grd.y, grd.v, grd.w);
        case 'decorrelation'
            [grd.v, grd.w] = l2grid(x_out2, y_out2, 1-cor_av, grd.x, grd.y, grd.v, grd.w);
        case 'intensity'
            [grd.v, grd.w] = l2grid(x_out1, y_out1, I_av1, grd.x, grd.y, grd.v, grd.w);
        case 'Chi_var'
            [grd.v, grd.w] = l2grid(x_out1, y_out1, X_var, grd.x, grd.y, grd.v, grd.w);
        case 'SI'
            [grd.v, grd.w] = l2grid(x_out1, y_out1, SI, grd.x, grd.y, grd.v, grd.w);
        case 'Sig_phi2'
            [grd.v, grd.w] = l2grid(x_out1, y_out1, sig_phi2, grd.x, grd.y, grd.v, grd.w);
        case 'Kp'
            [grd.v, grd.w] = l2grid(x_out1, y_out1, Kp, grd.x, grd.y, grd.v, grd.w);
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
    figure
    pcolorjw(xg,yg,covis.grid{2}.v);
    shading flat
    axis image;
    caxis([0 0.5]);
    colormap('jet');
    colorbar;
    xlabel('Easting of COVIS ( m )');
    ylabel('Northing of COVIS ( m )');
    hold on;
    plot(0,0,'.m','markersize',30);
    hold off;
    title('Decorrelation');


    figure
    pcolorjw(xg,yg,10*log10(covis.grid{3}.v));
    shading flat
    axis image;
    caxis([-50,-15]);
    colormap('jet');
    h = colorbar;
    title(h,'dB');
    xlabel('Easting of COVIS ( m )');
    ylabel('Northing of COVIS ( m )');
    hold on;
    plot(0,0,'.m','markersize',30);
    hold off;
    title('Target strength')

    figure
    pcolorjw(xg,yg,covis.grid{4}.v);
    shading flat
    axis image;
    caxis([0,2]);
    colormap('jet');
    colorbar;
    xlabel('Easting of COVIS ( m )');
    ylabel('Northing of COVIS ( m )');
    hold on;
    plot(0,0,'.m','markersize',30);
    hold off;
    title('log-amp fluctuation')

    figure
    si = covis.grid{5}.v;
    si(si==0) = nan;
    pcolorjw(xg,yg,si);
    shading flat
    axis image;
    caxis([0,1]);
    colormap('jet');
    colorbar;
    xlabel('Easting of COVIS ( m )');
    ylabel('Northing of COVIS ( m )');
    hold on;
    plot(0,0,'.m','markersize',30);
    hold off;
    title('Scintillation index');

    figure
    sp2 = covis.grid{6}.v;
    sp2(sp2==0) = nan;
    pcolorjw(xg,yg,sp2);
    shading flat
    axis image;
    caxis([0,2]);
    colormap('jet');
    colorbar;
    xlabel('Easting of COVIS ( m )');
    ylabel('Northing of COVIS ( m )');
    hold on;
    plot(0,0,'.m','markersize',30);
    hold off;
    title('Normalized Phase Variance');

    figure
    kp = covis.grid{7}.v;
    kp(kp==0) = nan;
    pcolorjw(xg,yg,kp);
    shading flat
    axis image;
    caxis([0,0.5]);
    colormap('jet');
    colorbar;
    xlabel('Easting of COVIS ( m )');
    ylabel('Northing of COVIS ( m )');
    hold on;
    plot(0,0,'.m','markersize',30);
    hold off;
    title('Normalized Amplitude Variance');
end
end
