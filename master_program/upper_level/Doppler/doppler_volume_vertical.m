% This program is used to calculate the vertical volume transport,
% e-folding radius, and centerline vertical flow rate of the plume

function [Q_out,Q_I,Q_d,U_out_uw,W,I_c,vr_back_m,vz_back_m,b_I,b_per,b_U,bw,center,vp] = doppler_volume_vertical(covis,zmin,zmax)
% Input:
% covis: output structure array from lower-level processing of the
%        Doppler-mode data
% zmin: lower bound of the vertical range ( m )
% zmax: upper bound of the vertical range ( m )
% Output:
% Q_out: Volume transport calculated using the method described in Xu et
%        al., 2013 ( m^3/s )
% Q_I: Volume transport calculated using the backscatter-intensity-weighted
%      method ( m^3/s )
% Q_d: Volume transport calculated using the method described in Xu et al., 2014 
% U_out_uw: centerline vertical flow rate calculated using the method
%           described in Xu et al., 2013 ( m/s )
% W: centerline vertical flow rate calculated using Darrell's method ( m/s )
% I_c: centerline backscatter intensity 
% vr_back_m: background radial velocity ( m/s )
% vz_back_m: backgorund vertical velocity ( m/s )
% b_I: e-folding radius calculated from horizontal backscatter intensity
%      cross-sections ( m )
% b_per: e-folding radius calculated from backscatter intensity
%      cross-sections perpendicular to the plume axis ( m )
% b_U: e-folding radius calculated from horizontal vertical velocity 
%      cross-sections ( m )
% bw: e-folding radius calculated using the method described in Xu et al.,
%     2014
% center: centerline coordinates ( m )
% vp: volume occupied by the plume within the selected vertical range ( m^3 ). 


flag =0;
I=covis.grid{2}.v;
I2 = covis.grid{2}.v_filt; % backscatter intensity cross-section perpendicular to plume axis
v=covis.grid{1}.v;
v_std = covis.grid{1}.std;
vr = covis.grid{1}.vr;
I_std = covis.grid{2}.std;
x = covis.grid{1}.x;
y = covis.grid{1}.y;
xx = squeeze(x(:,:,1));
yy = squeeze(y(:,:,1));



zs=10:0.5:60;
iz_min=find(zs==zmin);
iz_max=find(zs==zmax);

% locate the plume centerline
try
    [doppler_I]=doppler_intensity(covis,I,zmin,zmax,[-35,-20]);
catch me
    flag = 1;
    disp('an error happens while running doppler_intensity')
    Q_out = nan;
    U_out_uw = nan;
    U_out_w = nan;
    Q_I = nan;
    Q_d = nan;
    v_back_m = 0;
    b_I = 0;
    return
end
c_x = doppler_I.x0(:,2);
c_y = doppler_I.y0(:,2);
center(:,1) = c_x;
center(:,2) = c_y;
b_0 = doppler_I.radius(:,2);
I_0 = doppler_I.peak(:,2)+doppler_I.back(:,2);

try
    [doppler_std]=doppler_intensity(covis,I_std,zmin,zmax,[-35,-20]);
catch me
    disp('an error happens while running doppler_intensity')
end
try
    [doppler_I2]=doppler_intensity(covis,I2,zmin,zmax,[-35,-20]);
catch me
    doppler_I2.radius = nans(25,3);
    disp('an error happens while running dopppler_intensity')
end
% calculate the plume parameters and corresponding standard deviations
b_per=doppler_I2.radius(:,2); % radius of the perpendicular intensity profile (m)
b_I = doppler_I.radius(:,2); % radius of the horizontal intensity profile (m)
b_I_min=doppler_I2.radius(:,1); % lower boundary of the 95% confidence interval
b_I_max=doppler_I2.radius(:,3); % upper boundary of the 95% confidence interval
I_m = doppler_I.peak(:,2); % net centerline volume backscatter cross-section (m^-1)
I_b = doppler_I.back(:,2); % background volume backscatter cross-section
I_m_std = doppler_std.peak(:,2); % standard deviations for I_m
I_b_std = doppler_std.back(:,2); % stardard deivations for I_b
I_c = I_m+I_b; % centerline volume backscatter cross-section
I_c_std = sqrt(I_m_std.^2+I_b_std.^2); % standard deviation of the centerline volume backscatter cross-section
% Note that this relations is valid only if the
% background value is uncorrelated with the net
% value
% I_ave = (I_m.^2.*b_I.^2/2.*(1-exp(-2*R^2./b_I.^2))+I_b.^2*R^2-2*I_m.*I_b.*b_I.^2.*(exp(-R^2./b_I.^2)-1))...
%     ./(-b_I.^2.*I_m.*(exp(-R^2./b_I.^2)-1)+I_b.*R^2); % spatially averaged centerline back-scattering intensity
% I_ave_std = (I_m_std.^2.*b_I.^2/2.*(1-exp(-2*R^2./b_I.^2))+I_b_std.^2*R^2-2*I_m_std.*I_b_std.*b_I.^2.*(exp(-R^2./b_I.^2)-1))...
%     ./(-b_I.^2.*I_m_std.*(exp(-R^2./b_I.^2)-1)+I_b_std.*R^2); % spatially averaged intensity standard deviation

% % filter the data using one of the two filtering methods (region growing
% vs window function)
v_filt=zeros(size(v));
vr_back=zeros(size(v));
vz_back=zeros(size(v));
area = zeros(1,length(iz_min:iz_max));
%threshold = zeros(1,length(iz_min:iz_max));
index=0;
for nz=iz_min:iz_max
    index=index+1;
    [v_filt(:,:,nz),v_std(:,:,nz),area(index),vr_back(:,:,nz),vz_back(:,:,nz)]=doppler_filter(v(:,:,nz),v_std(:,:,nz),vr(:,:,nz),I(:,:,nz),doppler_I,index,'ccomp');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate average background velocity
vr_back_m = squeeze(nanmean(nanmean(vr_back,1),2));
vz_back_m = squeeze(nanmean(nanmean(vz_back,1),2));
vr_back_m = vr_back_m(iz_min:iz_max);
vz_back_m = vz_back_m(iz_min:iz_max);


%% Calculate volume flux using the method described in Xu et al. 2013
% spatial resolution of the grid
delta_x = 0.5;
delta_y = 0.5;
vv_filt = zeros(size(v_filt));
vv_std = zeros(size(v_filt));
vv_filt(:,:,iz_min:iz_max) = v_filt(:,:,iz_min:iz_max);
vv_std(:,:,iz_min:iz_max) = v_std(:,:,iz_min:iz_max);
% take out the velocity estimates smaller than 1 cm/s and larger than 0.8 m/s
i_v = find(abs(vv_filt(:))<1|abs(vv_filt(:))>80);
vv_filt(i_v)=0;
vv_std(i_v)=0;
% calculate the volume flux in the vertical direction by integrating the
% velocity estimates over horizontal cross-sections
Q_integral = 0.01*delta_x*delta_y*squeeze(nansum(nansum(vv_filt(:,:,iz_min:iz_max),2),1));
Q_integral_std = 0.01*delta_x*delta_y*sqrt(squeeze(nansum(nansum(vv_std(:,:,iz_min:iz_max).^2,2),1)));
Q_out(:,1) = Q_integral-Q_integral_std;
Q_out(:,2) = Q_integral;
Q_out(:,3) = Q_integral+Q_integral_std;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the volume flux using the backscatter cross-section integration method
% spatial resolution of the grid
delta_x = 0.5;
delta_y = 0.5;
lambda = 1.2; % concentration/velocity with ratio
I_o = I(:,:,iz_min:iz_max);
v_o = v(:,:,iz_min:iz_max);
%v_o_std = v_std(:,:,iz_min:iz_max);
% take out the velocity estimates smaller than 1 cm/s and larger than 0.8 m/s
i_v = find(abs(v_o(:))<1|abs(v_o(:))>80);
v_o(i_v)=0;
%v_o_std(i_v)=0;
Q_I = zeros(length(c_x),1);
for i = 1:length(c_x)
    r1 = sqrt((xx(:)-c_x(i)).^2+(yy(:)-c_y(i)).^2); % radial distance from the centerline
    b_w = b_0(i)/lambda;
    W1 = squeeze(v_o(:,:,i));
    I1 = squeeze(I_o(:,:,i));
    I01 = I_0(i);
    R = 4;
    ii1 = find(r1<R*b_w);
    W11 = W1(ii1);
    I11 = I1(ii1);
    WI_int = 0.01*delta_x*delta_y*sum(I11(:).*W11(:));
    Q_I(i) = WI_int*(1+lambda^2)/lambda^2/I01;
end



%% Calculate the volume flux,momentum flux using Darrell's method
thresh_Q = 0.1;   % Threshold for integral of velocity
thresh_P = thresh_Q^2;   % Threshold for integral velocity^2
windfact = 0.02;
e = 0;
Q_d = nan*ones(length(c_x),3);
P_d = nan*ones(length(c_x),3);
bw = nan*ones(length(c_x),3);
W = nan*ones(length(c_x),1);
vp = 0;
for nz = iz_min:iz_max
    e = e+1;
    bb = b_0(e);
    c_xx = c_x(e);
    c_yy = c_y(e);
    rr = sqrt((xx(:)-c_xx).^2+(yy(:)-c_yy).^2);
    ii = find(rr(:)<4*bb);
    if ~isempty(ii)
        v_sub = squeeze(v(:,:,nz));
        v_std_sub = squeeze(v_std(:,:,nz));
        I_sub = squeeze(I(:,:,nz));
        v_cut = v_sub(ii);
        I_cut = I_sub(ii);
        v_std_cut = v_std_sub(ii);
        x_cut = xx(ii);
        y_cut = yy(ii);
        v_cut = v_cut(:);
        I_cut = I_cut(:);
        x_cut = x_cut(:);
        y_cut = y_cut(:);
        X = [x_cut, y_cut];
        beta0_I=[max(I_cut(:)) bb*1e-6 0 c_xx*1e-6 c_yy*1e-6];
        try
            [beta_I,~,~] = nlinfit_xgy(X,I_cut,'gaussian2D',beta0_I);
            I_max = beta_I(1)+beta_I(3);
        catch
            I_max = max(I_cut(:));
        end
        % apply window filtering to velocity, based on intensity
        window = 1 - exp(-(I_cut/(windfact*I_max)).^4);
        vz_cut = v_cut.*window;
        vz_std_cut = v_std_cut.*window;
        beta0_W = [max(vz_cut(:)),bb*1e-6,0,c_xx*1e-6,c_yy*1e-6];
        try
            [beta_W,~,~] = nlinfit_xgy(X,v_cut,'gaussian2D',beta0_W);
            W_max = beta_W(1)+beta_W(3);
        catch
            W_max = max(vz_cut(:));
        end
        [Q_d(e,:),P_d(e,:),W(e)] = integ_gauss_guangyu(vz_cut,vz_std_cut,W_max,thresh_Q,thresh_P,0.5,0.5);
        int_ratio = Q_d(e,2).^2/P_d(e,2);
        Q_std = Q_d(e,2) - Q_d(e,1);
        P_std = P_d(e,2) - P_d(e,1);
        Q = Q_d(e,2);
        P = P_d(e,2);
        int_ratio_std = int_ratio*(2*Q_std/Q+P_std/P);
        bw2 = int_ratio/2/pi;
        bw(e,2) = sqrt(bw2); % mean radius
        bw(e,1) = bw(e,2)-1/2*bw(e,2)*int_ratio_std/int_ratio; % lower radius limit
        bw(e,3) = bw(e,2)+1/2*bw(e,2)*int_ratio_std/int_ratio; % upper radius limit
        % calculate the volume of the plume
        I_thresh = 10^-7; % backscatter cross-section threshold
        jj = find(I_cut>I_thresh);
        vp=vp+0.5^3*length(jj);
    else
        continue
    end
end

%% calculate the centerline vertical velocity and uncertainty

% unweighted Guaussian fit

try
    [doppler_W]=doppler_velocity(covis,v_filt,zmin,zmax,[]);
    %[doppler_W]=doppler_velocity(covis,v_filt,zmin,zmax,[-35,-20]);
catch me
    flag = 1;
    U_m_out = nan;
    disp('an error crops up while running doppler_velocity')
    return
end

%mean centerline velocity estimated by fitting a Gaussian curve to the
%velocity profiles perpendicular to the plume centerline.
U_m_out = (doppler_W.peak+doppler_W.back)/100;
U_m = U_m_out(:,2);
b_U = doppler_W.radius(:,2); % e-folding radius
U_ci_low = U_m_out(:,1); % lower boundary of the 95% CI
U_ci_high = U_m_out(:,3); % upper boundary of the 95% CI
x_c = doppler_W.x0(:,2); % centerline coordinate x
y_c = doppler_W.y0(:,2); % centerline coordiante y
j = 0;
U_m_std = zeros(1,length(iz_min:iz_max));
U_m_max = zeros(size(U_m_std));
U_m_min = zeros(size(U_m_std));
for nz = iz_min:iz_max
    j=j+1;
    v_std_sub = squeeze(v_std(:,:,nz));
    try
        %U_m_std(j)= interp2(covis_x(:,:,nz),covis_y(:,:,nz),v_std(:,:,nz),X(1,j),X(2,j))/100;
        F = TriScatteredInterp(xx(:),yy(:),v_std_sub(:))/100;
        U_m_std(j) = F(x_c(j),y_c(j));
    catch me
        U_m_std(j) = nan;
    end
    if ~isnan(U_ci_high(j))
        U_m_max(j) = U_m_std(j)+U_ci_high(j);
    else
        U_m_max(j) = U_m_std(j);
    end
    if ~isnan(U_ci_high(j))
        U_m_min(j) = -U_m_std(j)+U_ci_low(j);
    else
        U_m_min(j) = -U_m_std(j);
    end
end
% establish the output centerline velocity variable
U_out_uw = zeros(length(iz_min:iz_max),3);
U_out_uw(1:length(U_m_min),1)= U_m_min; %lower boundary
U_out_uw(1:length(U_m),2) = U_m; % mean value
U_out_uw(1:length(U_m_max),3) = U_m_max; % upper boundary

% weighted Guassian fit

% [W, X] = centerline_vel_b_nlin(covis,covis.grid{1}.v,20,32,[-35 -20], 'Gauss', 0);
% % use weighted nonlinear regression to calculate the centerline vertical
% % flow rate standard deviation
% [W_std] = centerline_vel_b_nlin(covis,covis.grid{1}.std,20,32,[-35,-20], 'Gauss', 0);
% W = W/100; % change the unit to m/s
% W_std = W_std/100;
% %U_m = W(2,:);
% %U_ci_high = W(1,:); % upper 95% CI boundary
% %U_ci_low = W(3,:); % lower 95% CI boundary
% %covis_x = covis.grid{1}.x;
% %covis_y = covis.grid{1}.y;
% %j=0;
% %[~,n] = size(X);
% %U_m_std = zeros(n,1);
% %U_m_max = zeros(n,1);
% %U_m_min = zeros(n,1);
% U_m_max = W+W_std;
% U_m_min = W-W_std;
% U_m = W;
%
%
% % establish the output centerline velocity variable
% U_out_w = zeros(length(iz_min:iz_max),3);
% U_out_w(1:length(U_m_min),1)= U_m_min; %lower boundary
% U_out_w(1:length(U_m),2) = U_m; % mean value
% U_out_w(1:length(U_m_max),3) = U_m_max; % upper boundary
end

% function for volume flux estimation using Darrell's method
function [Q_out, P_out, vc_old] = integ_gauss_guangyu(vz,vz_std,w_max,thresh_Q,thresh_P,delx,dely)
% Input:
% vz: vertical velocity estimates (cm/s)
% vz_std: standard deviation of velocity estimates (cm/s)
% thresh_Q: volume flux theshold
% thresh_P: momentum flux theshold
% delx: spatial interval in the x-axis direction
% dely: spatial interval in the y-axis direction
% Ouput:
% Q: volume flux (m^3/s);
% P: momentum flux
% vc_old: centerline vertical flow rate


% this function is adapted from Darrell's integ_gauss2.m by adding
% iteration scheme

% Does 2D integral of matrix_in assuming it is a 2D
% Gaussian (not necessarily isotropic), and only carrying
% the summation down to a specified level of the integrand
% thresh*max(max(matrix_in)). The tails of the integrand
% are accounted for by using

% int_-inf^inf int_-inf^inf exp(-x^2/bx^2 - y^2/by^2) dx dy =

%    1/(1-thresh)*int int  exp(-x^2/bx^2 - y^2/by^2) dx dy| exp(...) > thresh


%vc_old = max(vz(:)); % take the maximum velocity as the centerline velocity
vc_old = w_max;
error = 1;
iter = 0;
iter_max = 10;
while error>0.1&&iter<iter_max
    iter = iter +1;
    vz_norm = vz/vc_old;
    vz_std_norm = vz_std/vc_old;
    vz_keep = vz_norm(abs(vz_norm(:))>thresh_Q);
    vz_std_keep = vz_std_norm(abs(vz_norm(:))>thresh_Q);
    vz_std_keep_p = 2*vz_keep.*vz_std_keep;
    Q = 0.01*vc_old*delx*dely/(1-thresh_Q)*sum(vz_keep); % volume flux
    Q_std = 0.01*vc_old*delx*dely/(1-thresh_Q)*sqrt(sum(vz_std_keep(:).^2)); % volume flux uncertainty
    P = 0.01^2*vc_old^2*delx*dely/(1-thresh_P)*sum(vz_keep.^2); % specific momentum flux
    P_std = 0.01^2*vc_old^2*delx*dely/(1-thresh_P)*sqrt(sum(vz_std_keep_p(:).^2)); % momentum flux uncertainty
    %E = 1/2*0.01^3*vc_old^3*delx*dely/(1-thresh)*sum(vz_keep.^3); % energy flux
    vc_new = 2*P/Q*100;
    error = abs(vc_new - vc_old);
    vc_old = vc_new;
end
vc_old = vc_old/100;
Q_out = [Q-Q_std, Q, Q+Q_std];
P_out = [P-P_std, P, P+P_std];
if iter==iter_max
    disp('the solution failed to converge before the maximum number of interations is reached');
end
end