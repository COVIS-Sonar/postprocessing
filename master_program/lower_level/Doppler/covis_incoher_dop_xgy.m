function [vr_cov,vr_vel,rc,I_av,vr_std,I_std,covarsum] = covis_incoher_dop_xgy(hdr, dsp, range, data_burst)

% Coherent average of all input pings in burst is
% subtracted from all pings to reduce ground return through sidelobes.
%
% This function takes in beamformed data for one burst,
% and processes it into a 2D radial velocity slice,
% averaging over all pings in burst.
%
% output:
% vr_cov is the covariance ping-averaged line-of-sight velocity 
% vr_vel is the velocity ping-averaged line-of-sight velocity
% rc is the radial distance
% I_av is the ping-averaged volume backscattering coefficient 
% vr_std is the standard deviation of the radial velocity estimation.
% I_std is the standard deviation of the backscatter cross-section estimation.
% covarsum is the covariance function averaged over all the pings.
%---------------------------------------------------------------------
% versions - 
%  drj@apl.washington.edu
%  Edited 8/25/2010 by cjones@apl.washington.edu 
%  Edited 10/4/2010 by drj@apl
%  Edited 11/13/2019 by guangyux@uw.edu


cor = dsp.correlation;

% correlation range window size [number of samples]
if(~isfield(cor,'window_size'))
    cor.window_size = 0.001;
end
window_size = cor.window_size;

% correlation window overlap [number of samples]
if(~isfield(cor,'window_overlap'))
    cor.window_overlap = 0.4;
end
overlap = cor.window_overlap;

% ping-lag used to calculate covariance
lag = 8; 

sound_speed = hdr.sound_speed;
frequency = hdr.xmit_freq;
fsamp = hdr.sample_rate;
%pulse_width = hdr.pulse_width;

nwindow = round(window_size*fsamp);
noverlap = round(overlap*nwindow);
average = mean(data_burst,3);


% main loop
nbins=floor(size(data_burst,1)/(window-overlap));
covar = zeros(nbins,size(data_burst,2));
I1 = zeros(nbins,size(data_burst,2));
thetai = zeros(0);
for np = 1:size(data_burst,3)
    
    if(strfind(dsp.ping_combination.mode,'diff'))
        % Subtract average to reduce the contribution from the chimney
        sig_ping = data_burst(:,:,np) - average;
    else
        sig_ping = data_burst(:,:,np);
    end
    
    % Window in range
    [I1(:,:,np),~] = mag2_win(sig_ping,range,nwindow,noverlap);
    [covar(:,:,np),rc] = autocovar_win(sig_ping,range,nwindow,noverlap,lag);
    thetai(:,:,np)=angle(covar(:,:,np));
end    % End loop on np

I_av=sum(I1,3)/size(data_burst,3); % ping-averaged volume backscattering coefficient ( m^-1 )
covarsum=sum(covar,3);
thetac = angle(covarsum);
theta_std=std(thetai,0,3); % calculate the standard deviation of the angular frequency estimation
theta_m = mean(thetai,3); % calculate the mean of te angular frequency estimation
vr_vel = sound_speed*fsamp/(4*pi*frequency)*theta_m/lag; % velocity ping-averaged radial velocity ( m/s )
I_std=std(I1,0,3); % standard deviation of volume backscattering coefficient ( m^-1 )
vr_cov = sound_speed*fsamp/(4*pi*frequency)*thetac/lag; % covariance ping-averaged radial velocity ( m/s )
vr_std= sound_speed*fsamp/(4*pi*frequency)*theta_std/lag; % standard deviation of radial velocity ( m/s )
end






