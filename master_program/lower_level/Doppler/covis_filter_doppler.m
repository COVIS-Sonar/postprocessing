function [ data_out, filt, png ] = covis_filter_doppler( data_in, filt, png )
%
%  COVIS_FILER
%
% ----------
% Version 1.0 - cjones@apl.washington.edu 10/2010

% Set filter defaults
if(~isfield(filt,'bw'))
    filt.bw = 2;
end
if(~isfield(filt,'type'))
    filt.bw = 'butterworth';
end
if(~isfield(filt,'order'))
    filt.order = 4;
end
if(~isfield(filt,'decimation'))
    filt.decimation = 1;
end

% Apply Filter
if(strcmp(filt.status, 'on'))

    bw = filt.bw;     % Bandwdith is filt_bw/tau
    tau = png.hdr.pulse_width;     % Pulse length (sec)
    fsamp = png.hdr.sample_rate;   % Complex sampling frequency (Hz)
    order = filt.order;

    switch lower(filt.type)
        case {'butterworth'}
            [B,A] = butter(order, bw/fsamp/tau);
        otherwise
            disp('Unknown filter type');
    end
    
    % Lowpass filter
    %data_out = filter(B, A, data_in);
    data_out=filtfilt(B,A,data_in);
    % decimate
    R = filt.decimation;
    data_out = data_out(1:R:end,:);
    
    % update new sampling freq
    png.hdr.sample_rate = fsamp/R; 

end

end

