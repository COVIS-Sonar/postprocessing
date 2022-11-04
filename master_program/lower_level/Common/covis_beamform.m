function [bfm, bf_sig] = covis_beamform(bfm, raw_sig)
%
% covis_beamform.m
%
% Function to beamform selected Reson .dat file
%
% raw_sig =  matrix containing the baseband times series for each 
%            receiver element, should be of size [nsamp, nchans]            
% bfm.c = sound speed (m/s)
% bfm.fc = center frequency (Hz)
% bfm.fs = sampling frequency (Hz)
% bfm.first_samp = sample number of first data sample in raw_sig
% bfm.last_samp = sample number of last data sample in raw_sig
% bfm.type = type of beamforming ('fft', 'fast', or 'near')
% bfm.angle = beamforming angles (rad)
%
% Outputs:
% bf_sig =  matrix containing the baseband times series for each beam, 
%           size [nsamp, nbeams]
% bfm.range = range bins for each samples, size[nsamps,1] in meters
% bfm.angle = angle of each beam, size[1,nbeams] in radians
%
% -------------------------------------------------------------------------
% History
% version 1.0 - 31/1/10 DRJ@apl.washington.edu
% 05/10 cjones@apl.washington.edu - changed input/output arguments
% 06/10 cjones@apl.washington.edu - added fast beamformer
% 07/10 drj@apl.washington.edu - changed sign of phase in fast beamformer

global Verbose

if(~isfield(bfm,'fc'))
    fprintf('covis_beamform: error, fc must be defined\n');
    bf_sig = [];
    return;
end

% Set default type
if(~isfield(bfm,'type'))
   bfm.type = 'fast';
end

% Set default array length
if(~isfield(bfm,'array_length'))
   bfm.array_length = 0.408;  %  array length (m)
end

% Set default beams angles
if(~isfield(bfm,'start_angle'))
   bfm.start_angle = -64;
end
if(~isfield(bfm,'end_angle'))
   bfm.end_angle = 64;
end
if(~isfield(bfm,'num_beams'))
   if(bfm.fc == 200000) bfm.num_beams = 128; end
   if(bfm.fc == 396000) bfm.num_beams = 256; end
end

% Set default sound speed
if(~isfield(bfm,'c'))
   bfm.c = 1500;  % (m/s)
end

nsamps = size(raw_sig,1);
nchans = size(raw_sig,2);

% set default beam angles
% row vector containing nbeams angles at which beams are to be formed (deg)
nbeams = bfm.num_beams; % number of beams
start_angle = bfm.start_angle; 
end_angle = bfm.end_angle;
bfm.angle = (pi/180)*linspace(start_angle, end_angle, nbeams); % beam angles (rad)

% define beamformer params
type = bfm.type;
angle = bfm.angle;
L = bfm.array_length;  % array length (m)
c = bfm.c;  % sound speed (m/s)
fsamp = bfm.fs; % sampling freq (hz)
f = bfm.fc; % center frequency (hz)
lambda = c/f;  % wave length
k = 2*pi/lambda;  % wave number
first_samp = 1;

% slant range of each sample - round trip
dr = c/fsamp/2;
start_range = dr*first_samp;
end_range = dr*nsamps;
range = (start_range:dr:end_range)';

% set bfm structure range
bfm.range = range;

% w1 = column vector containing nch shading coefficients
w1 = hamming(nchans)';

lambda = c/f;  % wave length
omegac = 2*pi*f;
t_delay_max = L/2/c*sin(angle)';
delf = fsamp/nsamps;
%freqs =  delf*(-nsamps/2:nsamps/2-1);
%omegas = 2*pi*fftshift(freqs)';
freqs = delf*(0:nsamps-1);
omegas = 2*pi*freqs';
bf_sig = zeros(nsamps, nbeams);

if(Verbose > 2)
    fprintf('covis_beamform: type %s, f=%.2f, fs=%.2f, c=%.2f\n',type,f,fsamp,c);
    fprintf('covis_beamform: nsamps=%d, nbeams=%d\n',nsamps,nbeams);
    fprintf('covis_beamform: start_range=%.2f, end_range=%.2f\n',start_range, end_range);
end

% phase shift beamform
if (strcmp(type,'fast'))
    % phase shift to apply to each channel
    %  phi is a (nchans x nbeams) matrix
    d = (L/2)*linspace(-1,1,nchans)';  % element spacing
    phi = k*d*sin(angle); % phase shift
    % phase shift and sum to form beams by matrix multiplication
    bf_sig = raw_sig*(((w1')*ones(1,nbeams)).*exp(1i*phi));
    
    % fft time shift beamform
elseif (strcmp(type,'fft'))
    for nb = 1:nbeams
        t_ch = t_delay_max(nb)*linspace(-1,1,nchans);
        % Remove carrier phase shifts corresponding to assumed target azimuth
        %sf_phase_corr = fft((ones(nsamps,1)*exp(1i*omegac*t_ch)).*raw_sig.*(ones(nsamps,1)*w1));
        sf_phase_corr = fft(raw_sig); %xgy
        % Shift envelopes and sum over channels
        bf_sig(:,nb) = (sum((ifft(exp(1i*omegas*t_ch).*sf_phase_corr)).')).';
      
        %bf_sig(:,nb) = (sum((ifft(sf_phase_corr)).')).'; %xgy
    end
    
    % near field beamformer
elseif (strcmp(type,'near'))
    xs = linspace(-L/2,L/2,nchans);
    range = c/2/fsamp*[1:nsamps]';
    r_matrix = ranges*ones(1,nchans);
    x_matrix = ones(nsamps,1)*xs;
    for nb = 1:nbeams
        phi = pi*angles(nb)/180;
        t_ch = t_delay_max(nb)*linspace(-1,1,nchans);
        t_near = 1/c*(sqrt(r_matrix.^2+x_matrix.^2-2*sin(phi)*x_matrix.*r_matrix)...
            -r_matrix + c*ones(nsamps,1)*t_ch);
        % Make Fresnel phase corections for wavefront curvature
        s1_near = exp(1i*omegac*t_near).*raw_sig.*(ones(nsamps,1)*w1);
        % Remove carrier phase shifts for plane wave with assumed target azimuth
        sf_phase_corr = fft((ones(nsamps,1)*exp(1i*omegac*t_ch)).*s1_near);
        % Shift envelopes and sum over channels
        bf_sig(:,nb) = (sum((ifft(exp(1i*omegas*t_ch).*sf_phase_corr)).')).';
    end
    
else
    fprintf('covis_beamform: UNKNOWN BEAMFORMER TYPE/n');
end

if(Verbose > 4)
    fprintf('raw_sig min=%d, max=%d\n',min(min(abs(raw_sig).^2)),max(max(abs(raw_sig).^2)));
    fprintf('bf_sig min=%d, max=%d\n',min(min(abs(bf_sig).^2)),max(max(abs(bf_sig).^2)));
end


%figure(1); clf;
%subplot(2,1,1);
%plot(abs(bf_sig));
%subplot(2,1,2);
%plot(abs(raw_sig));
%pause;



