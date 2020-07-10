% This function is used to estimate the vertical velocity from the measured
% radial velocity component through the geometric conversion described in
% Jackson et al., 2003

function output = vertical_velocity(covis,z_fit_min,z_fit_max)
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

% set up grid parameters
% z_vent = 4.5; % height of Inferno vents above the base of COVIS ( m );
% z_fit_min = z_vent; % lower bound of vertical range of interest ( m )
% z_fit_max = 10; % upper bound of vertical range of interest ( m )
zg = squeeze(covis.grid{1}.z(1,1,:));
[~,iz1] = min(abs(zg-z_fit_min));
[~,iz2] = min(abs(zg-z_fit_max));
xg = squeeze(covis.grid{1}.x(1,:,1));
yg = squeeze(covis.grid{1}.y(:,1,1));

% calculate radial velocity from gridded covariance function
vr = covis.grid{1}.vr_cov;


% extract plume centerline properties
xs = -13.27; % x-coordinate of Inferno (m) (2019)
ys = -0.77; % y-coordinate of Inferno (m) (2019)
Ri = 4; % inner radius (m)
Ro = 8; % outer radius (m)
diagnosis = 0; % diagnosis toggle
grd = covis.grid{3};
doppler_I = get_plume_centerline_properties_fun(grd,xs,ys,Ri,Ro,z_fit_min,z_fit_max,diagnosis);

% coordiantes of the centerline
x0=doppler_I.x(:,2);
y0=doppler_I.y(:,2);
z0=doppler_I.z;
corr_ori(:,1)=x0;
corr_ori(:,2)=y0;
corr_ori(:,3)=z0;
x0_filt = movmean(x0,8,'omitnan');
y0_filt = movmean(y0,8,'omitnan');
x0_filt(isnan(x0)) = nan;
y0_filt(isnan(y0)) = nan;
corr_filt(:,1) = x0_filt;
corr_filt(:,2) = y0_filt;
corr_filt(:,3) = z0;

% Using a quadratic fit to parameterize the 3D centerline curve
t=(0:length(x0)-1)';
p1=polyfit(t,x0,2);
p2=polyfit(t,y0,2);
p3=polyfit(t,z0,2);
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
s=zeros(1,length(z0));
s(1)=0;
for i=2:length(z0)
s(i) = integral(F,t(1),t(i));
end



% calculate the linear-fit to the centerline (least-square);

[~,p,~] = best_fit_line(x0,y0,z0);
r=sqrt(p(1)^2+p(2)^2+p(3)^2);
inc_ang_lin=acos(abs(p(3))/r)*180/pi;
azimuth = atan2(p(2),p(1))*180/pi;

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
nping = [covis.burst(1:end).npings];
nping = median(nping); % number of pings in each burst
vr_std = covis.grid{1}.std/sqrt(nping);%standard deviation of the radial velocity averaged over npings
vr_var = vr_std.^2;
vz = zeros(size(vr));
vz_std = zeros(size(vr));
vz_var = zeros(size(vr));
x=covis.grid{1}.x;
y=covis.grid{1}.y;
z=covis.grid{1}.z;
vrc=zeros(length(z0),1); % centerline radial velocity
vrc_std=zeros(length(z0),1); % standard deviation of the centerline velocity
vc = zeros(length(z0),1);
v_h = zeros(length(z0),2);
v_h_std = zeros(length(z0),2);


% estimate the centerline radial velocity
j = 0;
for i = iz1:iz2
    j = j+1;
    vrc(j)=interp2(x(:,:,i),y(:,:,i),vr(:,:,i),x_fit(j),y_fit(j));
    vrc_std(j)=interp2(x(:,:,i),y(:,:,i),vr_std(:,:,i),x_fit(j),y_fit(j));
end

% smoothing the centerline axial velocity using Lowess smoothing
vrc_smooth = zeros(1,length(z0));
for j = 1:length(z0)
     z01 = z0(j);
     r = abs(z0 - z01);
     r = r/max(r);
     w = (1-r.^3).^3;
     E = [ones(length(z0),1), z0(:)-z01, (z0(:)-z01).^2];
     Ehat = diag(w)*E;
     vrc_sub_w = w(:).*vrc(:);
     a = Ehat\vrc_sub_w(:);
     vrc_smooth(j) = a(1); 
end

vrc = vrc_smooth;

j = 0;
for i=iz1:iz2
    j = j+1;
    % set a threshold on the angle between the plume centerline and the
    % acoustic line-of-sight to avoid negative centerline velocity
    if vrc(j)<0&&abs(ecerc(j))<0.2
        ecerc(j) = -0.2;
    elseif vrc(j)>=0&&abs(ecerc(j))<0.2
        ecerc(j) = 0.2;
    end
    % calculate the elevation angle sin_theta
    r=sqrt(x(:,:,j).^2+y(:,:,j).^2+z(:,:,j).^2);
    sin_theta=zg(i)./r;
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
Ig = covis.grid{2}.Id_filt;
cd_I = zeros(size(x));
alpha = zeros(1,length(z0));
j = 0;
for i=iz1:iz2
    j = j+1;
    [x_slice,y_slice]=meshgrid(xg,yg);
    hslice = surf(x_slice,y_slice,zeros(size(x_slice))+zg(i));
    %calculate rotation axis and angles
    theta = 90+180*atan2(dydt(j),dxdt(j))/pi;
    alpha(j) = 180*acos(dzdt(j))/pi;
    rotate(hslice,[theta,0],alpha(j),[corr_fit(j,1),corr_fit(j,2),corr_fit(j,3)]);
    %rotate(hslice,[theta,0],alpha(j));
    xd = get(hslice,'XData');
    yd = get(hslice,'YData');
    zd = get(hslice,'ZData');
    delete(hslice);
    h =figure(1);
    set(h,'Visible','off');
    hs1 = slice(x,y,z,10*log10(Ig),xd,yd,zd,'nearest');
    shading interp
    cd_I(:,:,i)=get(hs1,'CData'); 
    clear hs1
    cd_I(:,:,i)=10.^(cd_I(:,:,i)/10); % backscatter intensity over cross-sections perpendicular to the plume axis
end
cd_I(isnan(cd_I))=0;

% save results into a structure
output.iz = [iz1:iz2];
output.inc_ang_lin = inc_ang_lin;
output.azimuth = azimuth;
output.angles = angles;
output.ecerc = ecerc;
output.corr_ori = corr_ori;
output.corr_filt = corr_filt;
output.corr_fit = corr_fit;
output.I = grd.Id;
output.cd_I = cd_I;
output.v_h = v_h;
output.vc = vc;
output.vrc = vrc;
output.vz = vz;
output.vz_std = vz_std;
output.vr = vr;
output.s = s;
output.x = x;
output.y = y;
output.z = z;
end



