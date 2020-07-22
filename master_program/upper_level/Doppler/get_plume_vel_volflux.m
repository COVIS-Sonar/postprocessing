% This program is used to process gridded Doppler data to estimate plume
% vertical velocity and volume flux

%% calculate vertical velocities

z_min = 5;
z_max = 14;
z_lin_min = 5;
z_lin_max = 12;

% set up grid parameters
zg = squeeze(covis.grid{1}.z(1,1,:));
[~,iz1] = min(abs(zg-z_min));
[~,iz2] = min(abs(zg-z_max));
xg = squeeze(covis.grid{1}.x(1,:,1));
yg = squeeze(covis.grid{1}.y(:,1,1));



% extract plume centerline properties
xs = -13.27; % x-coordinate of Inferno (m) (2019)
ys = -0.77; % y-coordinate of Inferno (m) (2019)
Ri = 4; % inner radius (m)
Ro = 8; % outer radius (m)
diagnosis = 0; % diagnosis toggle
thresh = 3;
grd = covis.grid{3};
Ia = grd.Ia;
Id = grd.Id;
dI = Id./Ia;
Id(10*log10(dI)<thresh) = 10^-9;
grd.Id = Id;
doppler_I = get_plume_centerline_properties_fun(grd,xs,ys,Ri,Ro,z_min,z_max,diagnosis);

% coordiantes of the centerline
flag0 = doppler_I.Flag;
I0 = doppler_I.I(:,2);
x0=doppler_I.x(:,2);
y0=doppler_I.y(:,2);
z0=doppler_I.z;
corr_ori(:,1)=x0;
corr_ori(:,2)=y0;
corr_ori(:,3)=z0;
x02 = x0;
y02 = y0;
x02(flag0==1) = nan;
y02(flag0==1) = nan;
d = sqrt((x02-xs).^2+(y02-ys).^2);
x02(d>2*(z_max-z_min)) = nan;
y02(d>2*(z_max-z_min)) = nan;
flag0(isnan(x02)) = 1;
x0_filt = movmean(x02,8,'omitnan');
y0_filt = movmean(y02,8,'omitnan');
x0_filt(isnan(x02)) = nan;
y0_filt(isnan(y02)) = nan;
corr_filt(:,1) = x0_filt;
corr_filt(:,2) = y0_filt;
corr_filt(:,3) = z0;

% Using a quadratic fit to parameterize the 3D centerline curve
t=(0:length(x0_filt)-1)';
ii = find(~isnan(x0_filt));
p1=polyfit(t(ii),x0_filt(ii),2);
p2=polyfit(t(ii),y0_filt(ii),2);
p3=polyfit(t(ii),z0(ii),2);
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
covis_alt = 4.2;
ee=sqrt(x_fit.^2+y_fit.^2+(z_fit-covis_alt).^2);

% coordinates of the unit radial vector
xx=x_fit./ee;
yy=y_fit./ee;
zz=(z_fit-covis_alt)./ee;

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
ii = find(z0>z_lin_min&z0<z_lin_max&~isnan(x0_filt));
x0_sub = x0_filt(ii);
y0_sub = y0_filt(ii);
z0_sub = z0(ii);
[m,p,~] = best_fit_line(x0_sub-x0_sub(1),y0_sub-y0_sub(1),z0_sub-z0_sub(1));
t = -10:10;
x_lin_fit = m(1)+p(1)*t;
y_lin_fit = m(2)+p(2)*t;
z_lin_fit = m(3)+p(3)*t;
x_lin_fit = x_lin_fit+x0_sub(1);
y_lin_fit = y_lin_fit+y0_sub(1);
z_lin_fit = z_lin_fit+z0_sub(1);
corr_lin_fit(:,1)=x_lin_fit;
corr_lin_fit(:,2)=y_lin_fit;
corr_lin_fit(:,3)=z_lin_fit;
if p(3)<0
    p = -p;
end
r=sqrt(p(1)^2+p(2)^2+p(3)^2);
inc_ang_lin =acos(abs(p(3))/r)*180/pi;
azimuth_lin = atan2(p(2),p(1))*180/pi;

% calculate the inclination and declination angle based on quadratic fit.
inc_ang = acos(dzdt);
% calculate the azimuth angle based on the quadratic fit
dec_ang(dydt>=0) = acos(dxdt(dydt>=0)./sqrt(dxdt(dydt>=0).^2+dydt(dydt>=0).^2));
dec_ang(dydt<0) = -acos(dxdt(dydt<0)./sqrt(dxdt(dydt<0).^2+dydt(dydt<0).^2));

% calculate the angle between the sonar radial direction and the plume's
% axis
ang_sonar = acos(ecerc);

% assign the three different angles to angles
angles(:,1)=inc_ang*180/pi; % inclination
angles(:,2)=dec_ang*180/pi; % azimuth
angles(:,3)=ang_sonar*180/pi; % angle between sonar radial direction and plume axis

% Estimate the vertical velocity from the radial velocity
nping = [covis.burst(1:end).npings];
nping = median(nping); % number of pings in each burst
vr = covis.grid{1}.vr_cov(:,:,iz1:iz2);
vr_std = covis.grid{1}.std(:,:,iz1:iz2)/sqrt(nping);%standard deviation of the radial velocity averaged over npings
vr_var = vr_std.^2;
vz = zeros(size(vr));
vz_std = zeros(size(vr));
vz_var = zeros(size(vr));
xg=squeeze(covis.grid{1}.x(1,:,1));
yg=squeeze(covis.grid{1}.y(:,1,1));
zg=squeeze(covis.grid{1}.z(1,1,iz1:iz2));
vrc=zeros(length(z0),1); % centerline radial velocity
vrc_std=zeros(length(z0),1); % standard deviation of the centerline velocity
vc = zeros(length(z0),1);
v_h = zeros(length(z0),2);
v_h_std = zeros(length(z0),2);

% estimate the centerline radial velocity

for i = 1:length(z0)
    vrc(i)=interp2(xg,yg,vr(:,:,i),x_fit(i),y_fit(i));
    vrc_std(i)=interp2(xg,yg,vr_std(:,:,i),x_fit(i),y_fit(i));
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

[xxg,yyg] = meshgrid(xg,yg);
for j=1:length(z0)
    % set a threshold on the angle between the plume centerline and the
    % acoustic line-of-sight to avoid negative centerline velocity
    if vrc(j)<0&&abs(ecerc(j))<0.2
        ecerc(j) = -0.2;
    elseif vrc(j)>=0&&abs(ecerc(j))<0.2
        ecerc(j) = 0.2;
    end
    % calculate the elevation angle sin_theta
    r=sqrt(xxg.^2+yyg.^2+(zg(j)*ones(size(xxg))).^2);
    sin_theta=zg(j)./r;
    % estimate the vertical velocity component
    vz(:,:,j)=(vr(:,:,j)-vrc(j)*(dxdt(j)*xxg./r+dydt(j)*yyg./r)/ecerc(j))./sin_theta;
    C2 = sin_theta;
    A = 1./C2;
    vz_var(:,:,j) = A.^2.*vr_var(:,:,j);
    vz_std(:,:,j) = vz_var(:,:,j).^(1/2);
    % estimate the centerline velocity
    vc(j) = vrc(j)./ecerc(j);
    % estimate the horizontal velocity component
    v_h(j,1) = vrc(j)./ecerc(j)*dxdt(j); % eastward component
    v_h_std(j,1) = vrc_std(j)./ecerc(j)*dxdt(j); %standard deviation
    v_h(j,2) = vrc(j)./ecerc(j)*dydt(j); % northward component
    v_h_std(j,2) = vrc_std(j)./ecerc(j)*dydt(j); % standard deviation
end

% extract cross-section perpendicular to the plume's axis
Id_sub = Id(:,:,iz1:iz2);
cd_I = zeros(size(Id_sub));
alpha = zeros(1,length(z0));
for j=1:length(z0)
    [x_slice,y_slice]=meshgrid(xg,yg);
    hslice = surf(x_slice,y_slice,zeros(size(x_slice))+zg(j));
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
    hs1 = slice(xg,yg,zg,10*log10(Id_sub),xd,yd,zd,'nearest');
    shading interp
    cd_I(:,:,j)=get(hs1,'CData');
    clear hs1
    cd_I(:,:,j)=10.^(cd_I(:,:,j)/10); % backscatter intensity over cross-sections perpendicular to the plume axis
end
cd_I(isnan(cd_I))=10^-9;

%% calculate volume flux
dx = median(diff(xg));
dy = median(diff(yg));
windfact = 0.02;
Q = zeros(1,length(z0)); % volume flux ( m^3/s )
% filter velocities to eliminate noise outside the plume
R = 5;
area = nan(1,length(z0));
vz_filt = nan(size(vz));
for k=1:length(z0)
    if flag0(k) == 1
        continue
    end
    x01 = corr_fit(k,1);
    y01 = corr_fit(k,2);
    I01 = I0(k);
    I1 = Id_sub(:,:,k);

    vz1 = vz(:,:,k);
    ww = zeros(size(I1));
    ww(I1/I01>exp(-4)&abs(xxg-x01)<R&abs(yyg-y01)<R)=1;
    cc = bwconncomp(ww);
    object_prop = regionprops(cc);
    object_area = [object_prop.Area];
    [area_max,i_area]=max(object_area);
    ind = cc.PixelIdxList;
    ww_d = ww;
    ww_d(ind{i_area})=0;
    ww = ww-ww_d;
    area(k) = area_max*median(diff(xg))*median(diff(yg));
    window = 1 - exp(-(I1/(windfact*I01)).^4);
    vz_filt1 = vz1.*ww.*window;
    vz_filt(:,:,k) = vz_filt1;
    Q(k) = dx*dy*sum(vz_filt1(:))./mean(window(:));
end



% save results into a structure
output.iz = iz1:iz2;
output.inc_ang_lin = inc_ang_lin;
output.azimuth = azimuth_lin;
output.angles = angles;
output.ecerc = ecerc;
output.flag = flag0;
output.corr_ori = corr_ori;
output.corr_filt = corr_filt;
output.corr_fit = corr_fit;
output.corr_lin_fit = corr_lin_fit;
output.I = Id_sub;
output.Ip = cd_I;
output.v_h = v_h;
output.vc = vc;
output.vrc = vrc;
output.vz = vz;
output.vz_std = vz_std;
output.vr = vr;
output.Q = Q;
output.s = s;
output.x = xg;
output.y = yg;
output.z = zg;
%end



