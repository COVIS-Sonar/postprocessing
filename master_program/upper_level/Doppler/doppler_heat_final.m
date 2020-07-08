% This program is used to calculate the initial heat tranport of the plume
% above the North Tower of Grotto using the method described in Xu et al.,
% 2014

%% load the plume property time series and set up constant parameters

% directory and name of the plume property time series file
data_path = 'F:\COVIS\time_series\cat\final';
data_name = 'time_series_cat_default_correct_s.mat';
% load the time series 
load(fullfile(data_path,data_name));

% seawater constants obtained using the gsw toolbox based on the CTD 
% measurements of the plume from 2149 m to 2174 m depth
rho = 1.0373*10^3; % reference density (kg/m^3)
beta = 7.5010*10^-4;  % haline contraction coefficient (PSU^-1)
cp = 3.9200*10^3; % adiabetic heat capacity (J/kg/degree C)
alpha_T = 1.3180*10^-4; % thermal expansion coefficient ( degree C^-1);
g = 9.8; % graviational acceleration (m^2/s);
rho_cp_a_g = rho*cp/alpha_T/g; % W s^3/m^4
alpha_d = 0.095; % initial entrainment coefficient
depth_covis = 2193; % depth of the base of COVIS (m)

zs = 20:0.5:32;
z_cut = 24:0.5:32;
iz_min = find(zs == z_cut(1));
iz_max = find(zs == z_cut(end));

%% calculate heat flux by fitting the theoretical volume flux profile to the measured one
Q = squeeze(Q_d(iz_min:iz_max,2,:));
Q2 = squeeze(Q_d(:,2,:));
Q_g = squeeze(Q_out(:,2,:));
bU = squeeze(bw(iz_min:iz_max,2,:));
bU2 = squeeze(bw(:,2,:));
bI = b_I(iz_min:iz_max,:);
bI2 = b_I;
z_mid = mean(z_cut);
lambda = 1.06;
%lambda = 1;
Q(Q<0|Q>50) = nan;
Q2(Q2<0|Q2>50) = nan;
Q_g(Q_g<0|Q_g>50) = nan;
bU(bU<0|bU>20) = nan;
bU2(bU2<0|bU2>20) = nan;
bI(bI<0|bI>20) = nan;
bI2(bI2<0|bI2>20) = nan;
 ddate_t = [];
 HH_0_t = [];
 HH_0_e_t = [];
 zz_virtual_t = [];
 aalpha_t = [];
 aalpha_I_t = [];
 Aalpha_t = [];
 QQ_t = [];
 QQ2_t = [];
 QQ_g_t = [];
 bb_t = [];
 bb2_t = [];
 bbI_t = [];
 bbI2_t = [];
 ZZ_i_t = [];
 II_c_t = [];
 UU_t = [];
 UU_g_t = [];
 vp_t = [];
 error_Q = [];
 error_B = [];
       
[Y,M,D] = datevec(date);
e = 0;
for i_y = 2010:2015
    for i_m = 1:12
        jj = find(Y==i_y&M==i_m);
        if isempty(jj)
            continue;
        else
            e = e+1;
            date_seg = date(jj);
            Q_seg = Q(:,jj);
            Q2_seg = Q2(:,jj);
            Q_g_seg = Q_g(:,jj);
            bU_seg = bU(:,jj);
            bU2_seg = bU2(:,jj);
            bI_seg = bI(:,jj);
            bI2_seg = bI2(:,jj);
            I_c_seg = I_c(:,jj);
            %U_d_seg = squeeze(U_out(:,2,jj));
            U_d_seg = U_d(:,jj);
            vp_seg = vp(jj);
            U_g_seg = squeeze(U_out(:,2,jj));
            ang_inc_seg = ang_inc(jj);
            %ang_inc_seg = inc_ang_lin(jj);
            eta_U = zeros(1,length(date_seg));
            eta_I = zeros(1,length(date_seg));
            error_b = zeros(1,length(date_seg));
            zv_mid = zeros(1,length(date_seg));
            Z_i = zeros(length(z_cut),length(date_seg));
            z_virtual = zeros(1,length(date_seg));
            p2 = zeros(1,length(date_seg));
            Q_fit = zeros(length(z_cut),length(date_seg));
            H_0 = zeros(1,length(date_seg));
            error_q = zeros(1,length(date_seg));
            H_0_e = zeros(1,length(date_seg));
            p2_std = zeros(1,length(date_seg));
            alpha_dif = 1;
            ee = 0;
            while alpha_dif > 0.001
            for i = 1:length(date_seg)
                % perform linear regression on the plume radius profile to determine the
                % location of the virtual source
                bU_sub = bU_seg(:,i);
                bI_sub = bI_seg(:,i);
                Q_sub = Q_seg(:,i);
                p_U = polyfit(z_cut(~isnan(bU_sub))',bU_sub(~isnan(bU_sub)),1);
                p_I = polyfit(z_cut(~isnan(bI_sub))',bI_sub(~isnan(bI_sub)),1);
                eta_U(i) = p_U(1);
                eta_I(i) = p_I(1);
                b_fit_U = p_U(1)*z_cut(~isnan(bU_sub))+p_U(2);
                error_b(i) = sqrt(nanmean((b_fit_U(:)-bU_sub(~isnan(bU_sub))).^2))/nanmean(bU_sub);
                eta_c = alpha_d*6/5;
                alpha = alpha_d;
                zv_mid(i) = nanmean(bU_sub/eta_c);
                z_i = z_cut-z_mid+zv_mid(i);
                Z_i(:,i) = z_i;
                z_virtual(i) = depth_covis + p_U(2)/p_U(1);
                % linear fit
                E = z_i(~isnan(Q_sub)).^(5/3);
                p2(i) = E(:)\Q_sub(~isnan(Q_sub));
                if imag(p2(i))~=0
                   Q_fit(:,i) = nan;
                   H_0(i) = nan;
                   error_q(i) = nan;
                   H_0_e(i) = nan;
                else
                   Q_fit(:,i) = p2(i)*z_i.^(5/3);
                   b2 = (2*(5/(6*alpha))^4)/(3*pi^2*(1+lambda^2));
                   B_0 = p2(i)^3*b2; % initial buoyancy flux
                   H_0(i) = rho*cp/(g*alpha_T)*B_0/10^6; % initial heat flux
                   % calculate the uncertainty of the least square fit
                   n = Q_fit(:,i) - Q_sub; % fitting error
                   Error = sqrt(nanmean(n.^2))/nanmean(Q_seg(:,i)); % normalized RMS error
                   if isinf(Error)||imag(Error)~=0
                       error_q(i) = nan;
                   else
                       error_q(i) = Error;
                   end
                   sigma = std(n); % fitting error std
                   p2_std(i) = sqrt(sigma^2*inv(E*E')); % std of the slope
                   B_0_e = p2_std(i)^3*(2*(5/(6*alpha))^4)/(3*pi^2*(1+lambda^2)); % initial buoyancy flux uncertainty
                   H_0_e(i) = rho*cp/(g*alpha_T)*B_0_e/10^6; % initial heat flux uncertainty
                end
            end
            ii = find(ang_inc_seg*180/pi<10&error_q<0.2&error_b<0.05);
            if isempty(ii)
                alpha_d = 0.1083 + 0.01*randn(1,1);
                alpha_dif = 1;
                continue
            end
            ddate = date_seg(ii);
            HH_0 = H_0(ii);
            HH_0_e = H_0_e(ii);
            zz_virtual = z_virtual(ii);
            aalpha = eta_U(ii)/(6/5);
            aalpha_I = eta_I(ii)/(6/5);
            QQ = Q_seg(:,ii);
            QQ2 = Q2_seg(:,ii);
            QQ_g = Q_g_seg(:,ii);
            bb = bU_seg(:,ii);
            bb2 = bU2_seg(:,ii);
            bbI = bI_seg(:,ii);
            bbI2 = bI2_seg(:,ii);
            ZZ_i = Z_i(:,ii);
            II_c = I_c_seg(:,ii);
            UU = U_d_seg(:,ii);
            vvp = vp_seg(ii);
            UU_g = U_g_seg(:,ii);
            alpha_dif = abs(nanmean(aalpha) - alpha_d);
            alpha_d = nanmean(aalpha);
            
            end
            Aalpha = alpha_d*ones(1,length(aalpha));
        end
         ddate_t = cat(2,ddate_t,ddate); % time axis
         HH_0_t = cat(2,HH_0_t,HH_0); % heat transport time series ( MW )
         HH_0_e_t = cat(2,HH_0_e_t,HH_0_e); % error in heat transport estimate induced by nolinear fitting ( MW )
         zz_virtual_t = cat(2,zz_virtual_t,zz_virtual); % vertical source depth ( m )
         aalpha_t = cat(2,aalpha_t,aalpha); 
         aalpha_I_t = cat(2,aalpha_I_t,aalpha_I);
         Aalpha_t = cat(2,Aalpha_t,Aalpha);
         QQ_t = cat(2,QQ_t,QQ);
         QQ2_t = cat(2,QQ2_t,QQ2);
         QQ_g_t = cat(2,QQ_g_t,QQ_g);
         bb_t = cat(2,bb_t,bb);
         bb2_t = cat(2,bb2_t,bb2);
         bbI_t = cat(2,bbI_t,bbI);
         bbI2_t = cat(2,bbI2_t,bbI2);
         ZZ_i_t = cat(2,ZZ_i_t,ZZ_i);
         II_c_t = cat(2,II_c_t,II_c);
         UU_t = cat(2,UU_t,UU);
         vp_t = cat(2,vp_t,vvp);
         UU_g_t = cat(2,UU_g_t,UU_g);
         error_Q = cat(2,error_Q,error_q);
         error_B = cat(2,error_B,error_b);
    end
end

%% calculate the plume's volume based on its radius
[m,n] = size(QQ_t);
V = zeros(1,n);
for i = 1:m-1
    b1 = bb_t(i,:);
    b2 = bb_t(i+1,:);
    V = V+1/3*pi*(b1.^2+b1.*b2+b2.^2)*0.5;
end

%% calculate the virtual source depth based on the formula given in Hunt 2001
PP_t = QQ_t.^2./bb_t.^2/2/pi;% momentum flux
F0 = HH_0_t*10^6*g*alpha_T/cp/rho*2/pi;
M0 = PP_t(1,:)*2/pi;
Q0 = QQ_t(1,:)/pi;
M02 = PP_t*2/pi;
Q02 = QQ_t/pi;
F02 = repmat(F0,[17,1]);
Aalpha_t2 = repmat(Aalpha_t,[17,1]);
%T0 = HH_0_t./QQ_t(1,:)/cp/rho
L = 5*Q0.^2.*F0./4./Aalpha_t./M0.^(5/2); % source parameter
L2 = 5*Q02.^2.*F02./4./Aalpha_t2./M02.^(5/2);
phi = (L-1)./L;
% calculate the coefficients of phi
N = 5; % number of terms to use
delta = 0;
for j = 1:N;
coef = 3/5/(5^(j-1)*factorial(j)*(10*j-3))*prod(1+5*((1:j)-1));
delta = delta+coef*phi.^j;
end
z_avs = L.^(-1/5).*(1-delta);
Lm = 5./6./Aalpha_t.*(Q0./M0.^(1/2));
z_v = z_avs.*Lm; % depth of the virtual source below 24 m above the base of COVIS.
z_depth = depth_covis-24+z_v; % depth of the virtual source
ii = find(L<=0.5); % the condition underlying the formula is violated when L<=0.5
z_v(ii) = nan;
z_depth(ii) = nan;


% calculate N day average (without Oct 2010 data)
ii = find(ddate_t > datenum(2011,1,1,0,0,0));
date_sub = ddate_t(ii);
HH_sub = HH_0_t(ii);
HH_m = nanmean(HH_sub);
H_up = quantile(HH_sub,0.75);
H_low = quantile(HH_sub,0.25);
H_iq = H_up - H_low;
N = 0.4;
date_m = zeros(0);
H_m = zeros(0);
H_std = zeros(0);
date_cut = date_sub(1);
j = 0;
while date_cut < date_sub(end)
     j = j+1;
     ii = find(date_sub>=date_cut&date_sub<date_cut+N&~isnan(HH_sub));
     if length(ii)>=3
         date_m(j) = nanmean(date_sub(ii));
         H_m(j) = nanmean(HH_sub(ii));
         H_std(j) = nanstd(HH_sub(ii));
     else
         H_m(j) = nan;
         date_m(j) = nan;
         H_std(j) = nan;
     end
     date_cut = date_cut+N;
end
ii1 = find(~isnan(date_m)&~isnan(H_m)&H_std~=0);
date_m = date_m(ii1);
H_m = H_m(ii1);
H_std = H_std(ii1);


% calculate moving average of the time series
jj = find(ddate_t>datenum(2012,4,1,0,0,0)&ddate_t<datenum(2015,3,1,0,0,0));
%jj = find(ddate_t>datenum(2011,1,1,0,0,0)&ddate_t<datenum(2013,12,1,0,0,0));
[date_avg, H_avg, H_err,N] = heat_flux_movavg(ddate_t(jj), HH_0_t(jj), 1, 20, 30);



% fill the data gap in the heat flux time series with nans
[H_new,date_new] = gap_fill_nan_heat(HH_0_t,ddate_t,1/24/3600,3);


