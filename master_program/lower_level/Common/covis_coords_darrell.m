function [xw,yw,zw]=covis_coords_darrell(origin,range,azim,yaw,roll,pitch,central_head,declination)
% Input:
% origin: x-y-z coordinates of the sonar head
% range: sonar range ( m )
% azim: sonar azimuthal angles (radian)
% yaw: rotation relative to starting position in yaw (radian)
% roll: rotation relative to starting position in roll (radian)
% pitch: rotation relative to starting position in pitch (radian)
% central_head: central heading magnetic (degree)



% COVIS central pointing direction re true north
ang_rot = pi/180*(central_head + declination -360);

% rotation matrix for central heading
ROT = [cos(ang_rot) sin(ang_rot) 0 ; -sin(ang_rot) cos(ang_rot) 0 ; 0 0 1];


% Matrices giving coordinate transformations
R = [cos(roll) 0 -sin(roll) ; 0 1 0 ; sin(roll) 0 cos(roll)];
Y = [cos(yaw) -sin(yaw) 0 ; sin(yaw) cos(yaw) 0 ; 0 0 1];
P = [1 0 0 ; 0 cos(pitch) sin(pitch) ;  0 -sin(pitch) cos(pitch)];


% Receiver coordinates
xr = range(:)*sin(azim);
yr = range(:)*cos(azim);
zr = zeros(size(xr));
rr = [xr(:)';yr(:)';zr(:)'];

% World coordinates
rw = (R*P*Y)\rr;
rw = ROT*rw;
xw = reshape(rw(1,:),length(range),length(azim));
yw = reshape(rw(2,:),length(range),length(azim));
zw = reshape(rw(3,:),length(range),length(azim));


x0 = origin(1);
y0 = origin(2);
z0 = origin(3);

xw = xw+x0;
yw = yw+y0;
zw = zw+z0;



