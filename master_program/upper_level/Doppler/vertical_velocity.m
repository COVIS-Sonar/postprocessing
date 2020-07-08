% This function is used to estimate the vertical velocity from the measured
% radial velocity component through the geometric conversion described in
% Jackson et al., 2003

function [inc_ang_lin,azimuth,corr_ori,corr_fit,cd_I,v_h,vc,vz,vz_std,vr,vr_offset,ss] = vertical_velocity(covis,z_slice_min,z_slice_max,method)
% Intput:
% covis: output structure array from the lower-level processing of the
%        Doppler-mode data
% z_slice_min: lower bound of the vertical range ( m )
% z_slice_max: upper bound of the vertical range ( m )

% Output:
% inc_ang_lin: inclination angle of plume axis based on linear fit ( radian )
% azimuth: azimuth of plume axis based on linear fit ( radian )
% corr_ori: original plume center coordinates ( m )
% corr_ori: best-fit plume center coordinates ( m )
% cd_I: backscatter intensity over cross-sections perpendicular to the
%       plume axis
% v_h: horizontal cross-flow velocity ( cm/s )
% vc: centerline velocity along the tangential direction ( cm/s )
% vz: vertical velocity ( cm/s )
% vz_std: standard deviation of vertical velocity ( cm/s )
% vr: radial velocity with artificial offset removed ( cm/s )
% vr_offset: artifical radial velocity offset caused by clock drift ( cm/s )
% ss: distance from the vents along the plume axis ( m )

% Revome the artificial radial velocity offset caused by the clock drift
if strcmp(method,'covariance')
    disp('applying covariance gridding');
    [vr,vr_offset] = covis_covar_vr(covis);
    vr = vr - vr_offset; % radial velocity
else
    disp('applying velocity gridding');
    if isfield(covis.grid{1},'covar')
    [~,vr_offset] = covis_covar_vr(covis);
    vr=covis.grid{1}.vr-vr_offset; % radial velocity
    else
    vr = covis.grid{1}.vr;
    vr_offset = nan;
    end
end



% set the vertical range
zmin=10;
zmax=60;
space=covis.grid{2}.spacing;
range=covis.grid{2}.axis;
zs=range(5):space.dz:range(6);
iz_min=find(zs==zmin);
iz_max=find(zs==zmax);
iz_slice_min = find(zs == z_slice_min);
iz_slice_max = find(zs == z_slice_max);
Ig = covis.grid{2}.v; % volume backscatter cross-section from the output
nping = [covis.burst(1:end).npings];
nburst = length(nping); % number of bursts in the sweep
nping = median(nping); % number of pings in each burst
% locating the centerline of the plume
[doppler_I]=doppler_intensity(covis,Ig,zmin,zmax,[-35,-20]); 
% coordiantes of the centerline
x0=doppler_I.x0(:,2);
y0=doppler_I.y0(:,2);
z0=doppler_I.z;
corr_ori(:,1)=x0;
corr_ori(:,2)=y0;
corr_ori(:,3)=z0;


% Using a quadratic fit to parameterize the 3D centerline curve
t = 0:length(x0)-1;
z_fit_min = 20;
z_fit_max = 32;
p1=polyfit(t(z0>=z_fit_min&z0<=z_fit_max)',x0(z0>=z_fit_min&z0<=z_fit_max),2);
p2=polyfit(t(z0>=z_fit_min&z0<=z_fit_max)',y0(z0>=z_fit_min&z0<=z_fit_max),2);
p3=polyfit(t(z0>=z_fit_min&z0<=z_fit_max),z0(z0>=z_fit_min&z0<=z_fit_max),2);
% calculate the fitted quadratic curve
x_fit=p1(1)*t.^2+p1(2)*t+p1(3);
y_fit=p2(1)*t.^2+p2(2)*t+p2(3);
z_fit=p3(1)*t.^2+p3(2)*t+p3(3);
corr_fit(:,1)=x_fit;
corr_fit(:,2)=y_fit;
corr_fit(:,3)=z_fit;

% Using a third order polynomial to parameterize the 3-D plume centerline
% t = 0:length(x0)-1;
% z_fit_min = 20;
% z_fit_max = 32;
% iz = find(z0>=z_fit_min&z0<=z_fit_max);
% p1 = polyfit(t(iz)',x0(iz),3);
% p2 = polyfit(t(iz)',y0(iz),3);
% p3 = polyfit(t(iz),z0(iz),3);
% % calculate the third order polynomial fit
% x_fit=p1(1)*t.^3+p1(2)*t.^2+p1(3)*t+p1(4);
% y_fit=p2(1)*t.^3+p2(2)*t.^2+p2(3)*t+p2(4);
% z_fit=p3(1)*t.^3+p3(2)*t.^2+p3(3)*t+p3(4);
% corr_fit(:,1)=x_fit;
% corr_fit(:,2)=y_fit;
% corr_fit(:,3)=z_fit;

% caculate the unit radial vector of the fitted centerline
ee=sqrt(x_fit.^2+y_fit.^2+z_fit.^2);

% coordinates of the unit radial vector
xx=x_fit./ee;
yy=y_fit./ee;
zz=z_fit./ee;

% calculate the unit tangential vector of the fitted quadratic centerline
dxdt=2*p1(1)*t+p1(2);
dydt=2*p2(1)*t+p2(2);
dzdt=2*p3(1)*t+p3(2);
l=sqrt(dxdt.^2+dydt.^2+dzdt.^2);
%coordinates of the unit tangential vector
dxdt=dxdt./l;
dydt=dydt./l;
dzdt=dzdt./l;

% caculate the tangent vectors of the third order polynomial fit
% dxdt=3*p1(1)*t.^2+2*p1(2)*t+p1(3);
% dydt=3*p2(1)*t.^2+2*p2(2)*t+p2(3);
% dzdt=3*p3(1)*t.^2+2*p3(2)*t+p3(3);
% l=sqrt(dxdt.^2+dydt.^2+dzdt.^2);
% % coordinates of the unit tangential vector
% dxdt=dxdt./l;
% dydt=dydt./l;
% dzdt=dzdt./l;

% calculate the product of unit vectors of the plume's axis and acoustic
% line-of-sight
ecerc=dxdt.*xx+dydt.*yy+dzdt.*zz;



% calculate the fitted centerline distance
% quadratic fit
F = @(t)sqrt((2*p1(1)*t+p1(2)).^2+(2*p2(1)*t+p2(2)).^2+(2*p3(1)*t+p3(2)).^2); 
% cubic fit
%F =@(t)sqrt((3*p1(1)*t.^2+2*p1(2)*t+p1(3)).^2+(3*p2(1)*t.^2+2*p2(2)*t+p2(3)).^2+(3*p3(1)*t.^2+2*p3(2)*t+p3(3)).^2);
s1=zeros(1,length(z0));
s1(1)=0;
for i=2:length(z0)
s1(i) = quad(F,t(1),t(i));
end

% correct the centerline distance by subtracting the distance between COVIS
% and the top of the North Tower of Grotto, which is assumed to be 19 m
% above COVIS
s = s1 - s1(z0==19); % the top of the North Tower is at 17 m above the bottom of COVIS

% calculate the linear-fit to the centerline (least-square);
[~,vector,~]=best_fit_line(x0(z0>=z_fit_min&z0<=z_fit_max),y0(z0>=z_fit_min&z0<=z_fit_max),z0(z0>=z_fit_min&z0<=z_fit_max)');
if vector(3)<0;
    vector = vector*(-1); % reverse the vector if the z-component is negative
end
r=sqrt(vector(1)^2+vector(2)^2+vector(3)^2);
inc_ang_lin=acos(vector(3)/r); % calculate the inclination angle based on the linear fit;

% calculate the azimuth angle based on the linear fit
if(vector(1) < 0)
    azimuth=pi+atan(vector(2)/vector(1));
elseif(vector(1)>0 && vector(2)<0)
    azimuth=2*pi+atan(vector(2)/vector(1));
else
    azimuth= atan(vector(2)/vector(1)); 
end

% calculate the inclination and declination angle based on quadratic fit.
inc_ang = acos(dzdt);
% calculate the azimuth angle based on the quadratic fit
dec_ang(dydt>=0) = acos(dxdt(dydt>=0)./sqrt(dxdt(dydt>=0).^2+dydt(dydt>=0).^2));
dec_ang(dydt<0) = -acos(dxdt(dydt<0)./sqrt(dxdt(dydt<0).^2+dydt(dydt<0).^2));

% calculate the angle between the sonar radial direction and the plume's
% axis
ang_sonar = acos(ecerc);

% assign the three different angles to angles
angles(:,1)=inc_ang; % inclination
angles(:,2)=dec_ang; % azimuth
angles(:,3)=ang_sonar; % angle between sonar radial direction and plume axis

% Estimate the vertical velocity from the radial velocity
vr_std = covis.grid{1}.std/sqrt(nping);%standard deviation of the radial velocity averaged over npings
vr_var = vr_std.^2;
vz = zeros(size(vr));
vz_std = zeros(size(vr));
vz_var = zeros(size(vr));
x=covis.grid{1}.x;
y=covis.grid{1}.y;
z=covis.grid{1}.z;
vrc=zeros(1,length(x0)); % centerline radial velocity
vrc_std=zeros(1,length(x0)); % standard deviation of the centerline velocity
v_h = zeros(length(iz_min:iz_max),2);
v_h_std = zeros(size(v_h));
vc = zeros(vrc);
j=0;

% estimate the centerline radial velocity
for i = iz_min:iz_max
    j=j+1;
    vrc(j)=interp2(x(:,:,i),y(:,:,i),vr(:,:,i),x_fit(j),y_fit(j));
    vrc_std(j)=interp2(x(:,:,i),y(:,:,i),vr_std(:,:,i),x0(j),y0(j));
end

% smoothing the centerline axial velocity using Lowess smoothing
z_cut = z_slice_min:0.5:z_slice_max;
vrc_sub = vrc(iz_slice_min:iz_slice_max);
vrc_smooth = zeros(1,length(z_cut));
for i = 1:length(z_cut)
     z0 = z_cut(i);
     r = abs(z_cut - z0);
     r = r/max(r);
     w = (1-r.^3).^3;
     E = [ones(length(z_cut),1), z_cut(:)-z0, (z_cut(:)-z0).^2];
     Ehat = diag(w)*E;
     vrc_sub_w = w.*vrc_sub;
     a = Ehat\vrc_sub_w(:);
     vrc_smooth(i) = a(1); 
end

vrc(iz_slice_min:iz_slice_max) = vrc_smooth;
j=0;
for i=iz_min:iz_max
    j=j+1;
    % set a threshold on the angle between the plume centerline and the
    % acoustic line-of-sight to avoid negative centerline velocity
    if vrc(j)<0&&abs(ecerc(j))<0.2
        ecerc(j) = -0.2;
    elseif vrc(j)>=0&&abs(ecerc(j))<0.2
        ecerc(j) = 0.2;
    end
    % calculate the elevation angle sin_theta
    r=sqrt(x(:,:,i).^2+y(:,:,i).^2+z(:,:,i).^2);
    sin_theta=zs(i)./r;
    % estimate the vertical velocity component
    vz(:,:,i)=(vr(:,:,i)-vrc(j)*(dxdt(j)*x(:,:,i)./r+dydt(j)*y(:,:,i)./r)/ecerc(j))./sin_theta;
    C2 = sin_theta;
    A = 1./C2;
    vz_var(:,:,i) = A.^2.*vr_var(:,:,i);
    vz_std(:,:,i) = vz_var(:,:,i).^(1/2);
    % estimate the centerline velocity
    vc(j) = vrc(j)./ecerc(j);
    % estimate the horizontal velocity component
    v_h(j,1) = vrc(j)./ecerc(j)*dxdt(j); % eastward component
    v_h_std(j,1) = vrc_std(j)./ecerc(j)*dxdt(j); %standard deviation
    
    v_h(j,2) = vrc(j)./ecerc(j)*dydt(j); % northward component
    v_h_std(j,2) = vrc_std(j)./ecerc(j)*dydt(j); % standard deviation
    
end

% extract cross-section perpendicular to the plume's axis
xx=-40:0.5:-5;
yy=-25:0.5:5;
cd_I = zeros(61,71,101);
alpha = zeros(size(iz_slice_min:1:iz_slice_max));
for i=iz_slice_min:1:iz_slice_max
 
    [x_slice,y_slice]=meshgrid(xx,yy);
    hslice = surf(x_slice,y_slice,zeros(size(x_slice))+zs(i));
    %calculate rotation axis and angles
    theta = 90+180*atan2(dydt(i),dxdt(i))/pi;
    alpha(i) = 180*acos(dzdt(i))/pi;
    rotate(hslice,[theta,0],alpha(i),[corr_fit(i,1),corr_fit(i,2),corr_fit(i,3)]);
    %rotate(hslice,[theta,0],alpha(j));
    xd = get(hslice,'XData');
    yd = get(hslice,'YData');
    zd = get(hslice,'ZData');
    delete(hslice);
    h =figure(1);
    set(h,'Visible','off');
    hs1 = slice(covis.grid{2}.x,covis.grid{2}.y,covis.grid{2}.z,10*log10(Ig),xd,yd,zd,'nearest');
    shading interp
    cd_I(:,:,i)=get(hs1,'CData'); 
    clear hs1
    cd_I(:,:,i)=10.^(cd_I(:,:,i)/10); % backscatter intensity over cross-sections perpendicular to the plume axis
end
cd_I(isnan(cd_I))=0;
ss = s(iz_slice_min:1:iz_slice_max); % axial distance between vertical height z_slice_min and z_slice_max
end


%% This function is used to calculate the radial velocity and radial 
% velocity offset from the covariance functions. Note that the offset is 
% caused by the timing mismatch between the clocks in COVIS for signal
% generation and digitization.

function [vr_covar,vr_offset] = covis_covar_vr(covis)
% input:
% covis: covis structure array.
% output:
% va_covar: radial velocity (cm/s)
% vr_offset: offset in the radial velocity estimates (cm/s).
covar = covis.grid{1}.covar;
covar_offset = covis.grid{1}.offset_covar;
fs = covis.processing.beamformer.fs; % sampling rate (Hz);
c = covis.processing.beamformer.c; % sound speed (m/s);
fc = covis.processing.beamformer.fc; % central frequency (Hz);
thetac = angle(covar_offset);
thetaci = angle(covar);
vr_offset = 100*c*fs/(4*pi*fc)*thetac/8;
vr_covar = 100*c*fs/(4*pi*fc)*thetaci;
end

