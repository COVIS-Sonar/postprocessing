function  data2_corr = covis_phase_correct(ping, data1, data2)
%
% Corrects phase jitter of data2 by comparing the transmitted 
% signals recorded in a selected channel (chn) via cross-talk.  
 
% originally written by drj@apl.washington.edu 29 July 2010
% edited by guangyux@uw.edu in Aug 2019 to use the source waveform recorded
% in a regular channel through cross-talk for phase correction. This change
% was made in response to the malfunction of the monitor channel. 

chn = size(data1,2); % # of the reference channel


sample_rate = ping.hdr.sample_rate;
pulse_width = ping.hdr.pulse_width;

% Number of samples to keep from monitor channel
nkeep = round(3*sample_rate*pulse_width);

% The monitor channel is the last channel
monitor1 = data1(1:nkeep,chn);
monitor2 = data2(1:nkeep,chn);

% Cross-correlate monitor signals to determine phase jitter
phase_jitter = angle(sum(monitor2.*conj(monitor1)));

% Correct phase of data2 to agree with data1
data2_corr = exp(-1i*phase_jitter)*data2;

