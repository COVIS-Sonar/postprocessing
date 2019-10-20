function [xw,yw,zw]=covis_coords(origin,rs,azim,compass,roll,elev)
%
% Transform covis sonar coords (origin,rs,azim) into rectangular world coords (xw,yw,zw) 
% given the sonar attitude defined by the TCM sensor (compass,roll,elevation).  
% origin = [x0,y0,z0], rs = slant_range, azim = beam_angle.
%
% The output world coord (xw,yw,zw) are defined as:
%   X-direction is East.
%   Y-direction is North.
%   Z-direction is Up from the seafloor.
%
% All angles are given in radians.
%
% Inputs:
%  rs must be a column vector or a scalar,
%  azim must be a row vector or a scalar, 
%  compass, roll, and elev angles must be scalar.
%
% Returns a matrix for each of the sonar rectangular coords xw, yw, and zw.
% Each matrix is of size [length(rs),length(azim)], corresponding to the 
% x, y, and z coordinates of the each input range (rs) and azimuth (azim). 
%
% ---
% Version 1.0 - cjones@apl.washington.edu 08/2010
%



x0 = origin(1);
y0 = origin(2);
z0 = origin(3);

% transform from sensor [r,azim] to sensor [xs,ys,zs]
xs = rs*sin(azim);
ys = rs*cos(azim)*cos(elev);
zs = rs*cos(azim)*sin(elev);

% Roll (rotation around ys axis)
xs =  xs*cos(roll) + zs*sin(roll);
zs = -xs*sin(roll) + zs*cos(roll);

% compass and translation
xw = x0 + xs*cos(compass) + ys*sin(compass);
yw = y0 - xs*sin(compass) + ys*cos(compass);
zw = z0 + zs;

return;


% % The following is for Euler angles
% 
% % change angle signs to match the sonar
% roll = -roll;
% pitch = pitch;
% yaw = yaw;
% 
% % reshape the coords
% Rs=[xs(:)'; ys(:)'; zs(:)'];
% 
% % transform into world coords [xw,yw,zw]
% C1=[ cos(yaw) sin(yaw) 0; 
%     -sin(yaw) cos(yaw) 0; 
%      0        0        1];
% C2=[ 1  0         0;
%      0  cos(roll) sin(roll); 
%      0 -sin(roll) cos(roll)];
% C3=[ cos(pitch) 0 -sin(pitch); 
%      0          1  0;
%      sin(pitch) 0  cos(pitch)];
% 
% C = C1*C2*C3;  % pitch*roll*yaw transform order
% R = C*Rs; % do the transform
% 
% % reshape into [length(rs),length(azim)]
% xw = x0 + reshape(R(1,:),size(xs));
% yw = y0 + reshape(R(2,:),size(ys));
% zw = z0 + reshape(R(3,:),size(zs));
% 
% 
