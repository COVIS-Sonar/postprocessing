% This program is used to concatanate the monthly plume property time
% series into a continuous multi-year-long time series

% directory of the monthly time series
path_to_data = 'F:\COVIS\time_series\final_time_series\default_correct_s_final';
seq_file = dir(fullfile(path_to_data,'*.mat'));
file1 = load(fullfile(path_to_data,seq_file(1).name),'Q_out','Q_I','Q_d','I_c','U_out_uw','U_d','v_back','vr_back','vr_offset','inc_ang_lin','azimuth','date','s','b_I','bw','b_U','vp');
date = file1.date;
Q_out = file1.Q_out;
U_out = file1.U_out_uw;
Q_I = file1.Q_I;
Q_d = file1.Q_d;
U_d = file1.U_d;
ang_inc = file1.inc_ang_lin;
ang_azi = file1.azimuth;
vp = file1.vp;
s = file1.s;
b_I = file1.b_I;
bw = file1.bw;
b_U = file1.b_U;
if isfield(file1,'v_back')
   v_back = file1.v_back;
else
    v_back = file1.vr_back;
end
vr_offset = file1.vr_offset;
I_c = file1.I_c;

for i = 2:length(seq_file)
    file = load(fullfile(path_to_data,seq_file(i).name),'Q_out','Q_I','Q_d','U_out_uw','U_d','I_c','v_back','vr_back','vr_offset','inc_ang_lin','azimuth','date','s','b_I','bw','b_U','vp');
    Q_out = cat(3,Q_out,file.Q_out);
    Q_I = cat(2,Q_I,file.Q_I);
    Q_d = cat(3,Q_d,file.Q_d);
    U_out = cat(3,U_out,file.U_out_uw);
    U_d = cat(2,U_d,file.U_d);
 
    if isfield(file,'v_back')
       v_back = cat(2,v_back,file.v_back);
    else
       v_back = cat(2,v_back,file.vr_back);
    end
    vr_offset = cat(2,vr_offset,file.vr_offset);
    I_c = cat(2,I_c,file.I_c);
    ang_inc = cat(2,ang_inc,file.inc_ang_lin);
    ang_azi = cat(2,ang_azi,file.azimuth);
    date = cat(2,date,file.date);
    s = cat(2,s,file.s);
    b_I = cat(2,b_I,file.b_I);
    bw = cat(3,bw,file.bw);
    b_U = cat(2,b_U,file.b_U);
    vp = cat(2,vp,file.vp);
end
[date,I]=sort(date);
Q_out = Q_out(:,:,I);
Q_d = Q_d(:,:,I);
Q_I = Q_I(:,I);
U_out = U_out(:,:,I);
U_d = U_d(:,I);
v_back = v_back(:,I);
vr_offset = vr_offset(I);
I_c = I_c(:,I);
ang_inc = ang_inc(I);
ang_azi = ang_azi(I);
s = s(I);
vp = vp(I);
b_I = b_I(:,I);
bw = bw(:,:,I);
b_U = b_U(:,I);

% eliminate the redundant samples
ii = find(diff(date)==0);
if ~isempty(ii)
ee = [1:ii(1)-1,ii,ii(end)+2:length(date)];
date = date(ee);
Q_out = Q_out(:,:,ee);
Q_I = Q_I(:,ee);
Q_d = Q_d(:,:,ee);
U_out = U_out(:,:,ee);
U_d = U_d(:,ee);
ang_inc = ang_inc(ee);
ang_azi = ang_azi(ee);
s = s(ee);
vp = vp(ee);
b_I = b_I(:,ee);
bw = bw(:,:,ee);
b_U = b_U(:,ee);
end

%% calculate the vertically averaged volume flux and centerline vertical
% flow rate
U = squeeze(U_out(:,2,:));
Q = squeeze(Q_out(:,2,:));
U_m_v = nanmedian(U,1);
Q_m_v = nanmedian(Q,1);
ii_u = find(~isnan(U_m_v));
ii_q = find(~isnan(Q_m_v));
date_u = date(ii_u);
date_q = date(ii_q);
U_m_v = U_m_v(ii_u);
Q_m_v = Q_m_v(ii_q);

%% perform N-day average over the vertically averaged volume flux and
% centerline vertical flow rate
N = 30;

% centerline vertical flow rate
date_m_n_u = zeros(0);
U_m_n = zeros(0);
date_cut_u = date_u(1);
j = 0;
while date_cut_u < date_u(end)
     j = j+1;
     ii = find(date_u>=date_cut_u&date_u<date_cut_u+N);
     if ~isempty(ii)
         date_m_n_u(j) = mean(date_u(ii));
         U_m_n(j) = mean(U_m_v(ii));
     else
         U_m_n(j) = nan;
         date_m_n_u(j) = nan;
     end
     date_cut_u = date_cut_u+N;
end

% volume flux
date_m_n_q = zeros(0);
Q_m_n = zeros(0);
date_cut_q = date_q(1);
j = 0;
while date_cut_q < date_q(end)
     j = j+1;
     ii = find(date_q>=date_cut_q&date_q<date_cut_q+7);
     if ~isempty(ii)
         date_m_n_q(j) = mean(date_q(ii));
         Q_m_n(j) = mean(Q_m_v(ii));
     else
         Q_m_n(j) = nan;
         date_m_n_q(j) = nan;
     end
     date_cut_q = date_cut_q+N;
end

%% perform monthly average on the vertically averaged centerline vertical flow rate and volume flux time series

% centerline vertical flow rate
V_u = datevec(date_u);
year = sort(unique(V_u(:,1)));
month = sort(unique(V_u(:,2)));
U_m_m = zeros(0);
Q_m_m = zeros(0);
date_m_m = zeros(0);
e = 0;
for i = 1:length(year)
    for j = 1:length(month)
        d_low = datenum(year(i),month(j),1,0,0,0);
        d_up = datenum(year(i),month(j),31,23,59,59);
        ii_u = find(date_u>=d_low&date_u<=d_up);
        ii_q = find(date_q>=d_low&date_q<=d_up);
        if ~isempty(ii_u)
            e = e+1;
            U_m_m(e) = nanmean(U_m_v(ii_u));
            Q_m_m(e) = nanmean(Q_m_v(ii_q));
            date_m_m(e) = d_low;
        else
            continue;
        end
    end
end