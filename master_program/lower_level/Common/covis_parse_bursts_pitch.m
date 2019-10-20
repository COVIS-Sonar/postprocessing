%  This function searches through a list of pings and find the 
%  consequtive pings that have the same elevation value.  Each set of
%  consequtive pings in called a burst.
function [ burst ] = covis_parse_bursts_pitch(png)
% Input:   
% png: a structure array that includes sonar orientation angles and
% corresponding motor positions
% Output:
% burst: a structure array with the following fields:
%   burst.elev = elevation angle of burst
%   burst.npings = number of pings in burst
%   burst.start_ping = start ping number in burst
%
% ----------
% Version 1.0 - guangyux@uw.edu (based on the original code written by
% Chris Jones in 2010)

global Verbose;
pitch_rot = [png.rot_pitch];
pitch_tcm = [png.sen_pitch];
ping_no = [png.num];
[~,ia,~] = unique(pitch_rot,'stable');
for i = 1:length(ia)
    burst(i).pitch = pitch_tcm(ia(i));
    if i ~= length(ia)
        burst(i).ping = ping_no(ia(i):ia(i+1)-1);
        burst(i).npings = ia(i+1)-ia(i);
    else
        burst(i).ping = ping_no(ia(i):length(ping_no));
        burst(i).npings = length(ping_no)-ia(i)+1;
    end
    if(Verbose > 1)
        fprintf('Burst %d: pitch=%f, npings=%d, start_ping=%d\n', ...
            i, burst(i).pitch, burst(i).npings, burst(i).ping(1))
    end
end


end

