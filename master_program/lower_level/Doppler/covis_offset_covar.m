% This function is used to calculate the covariance function of the monitor
% signal. The covariance function will be used to correct the offset in the
% radial velocity estimates caused by the timing mismatch between the
% clocks in COVIS for data generation and digitization.
function [covar] = covis_offset_covar(data,ping)
% input:
% data: baseband complex signal
% ping: structure array of the constant parameters for each ping
% output:
% covar: covariance function
nchan = size(data,2);

sample_rate = ping.hdr.sample_rate;
pulse_width = ping.hdr.pulse_width;


% Number of samples to keep from monitor channel
nkeep = round(3*sample_rate*pulse_width);

% The monitor channel is the last channel
monitor = data(1:nkeep,nchan);

% calculate the covariance function
covar = sum(monitor(1:end-8).*conj(monitor(9:end)));




end