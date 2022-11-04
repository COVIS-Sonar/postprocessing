function [vr_cov,vr_vel,vr_std,I_a,I_d,covarsum,rc] = covis_incoher_dop_xgy(hdr, dsp, range, data_burst)

cor = dsp.correlation;


% correlation range window size [number of samples]

window_size = cor.window_size;

% correlation window overlap [number of samples]
overlap = cor.window_overlap;

% ping-lag used to calculate covariance
if(~isfield(cor,'nlag'))
    cor.nlag = 8;
end
lag = cor.nlag;

sound_speed = hdr.sound_speed;
frequency = hdr.xmit_freq;
fsamp = hdr.sample_rate;
pulse_width = hdr.pulse_width;

nwindow = round(window_size*fsamp);
noverlap = round(overlap*nwindow);
average = mean(data_burst,3);


% main loop
nbins=floor(size(data_burst,1)/(nwindow-noverlap));
covar = zeros(nbins,size(data_burst,2),size(data_burst,3));
I1_d = zeros(size(covar));
thetai = zeros(size(covar));
for np = 1:size(data_burst,3)
    sig_ping_d = data_burst(:,:,np) - average;
    [I1_d(:,:,np),~] = mag2_win(sig_ping_d,range,nwindow,noverlap);
    [covar(:,:,np),rc] = autocovar_win(data_burst(:,:,np),range,nwindow,noverlap,lag);
    thetai(:,:,np)=angle(covar(:,:,np));
end    % End loop on np
[I_a,~] = mag2_win(average,range,nwindow,noverlap); % ping-averaged volume backscattering coefficient ( m^-1 )
I_d= mean(I1_d,3); % ping-differenced volume backscattering coefficient ( m^-1 )
covarsum=nansum(covar,3);
thetac = angle(covarsum);
theta_std=nanstd(thetai,0,3); % calculate the standard deviation of the angular frequency estimation
theta_m = nanmean(thetai,3); % calculate the mean of te angular frequency estimation
vr_vel = sound_speed*fsamp/(4*pi*frequency)*theta_m/lag; % velocity ping-averaged radial velocity ( m/s )
vr_cov = sound_speed*fsamp/(4*pi*frequency)*thetac/lag; % covariance ping-averaged radial velocity ( m/s )
vr_std= sound_speed*fsamp/(4*pi*frequency)*theta_std/lag; % standard deviation of radial velocity ( m/s )
end






