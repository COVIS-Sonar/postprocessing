% This function is used to build the time series of plume properties 
% during the selected time range
function covis_time_series(year,month,dates,z_min,z_max)
% Inputs:
% year: (e.g., 2014)
% month:(e.g., May)
% dates: a vector containing the days in a month (e.g., 1:30 for the 1st to
%        the 30th)
% z_min: lower bound of the vertical range ( m )
% z_max: upper bound of the vertical range ( m )

% directory of the Doppler data
file_path = ['F:\COVIS\covis_doppler_data\',year,'\'];
zs = 10:0.5:60;
i_z = find(zs>=z_min&zs<=z_max);
% directory under which the combined and geometrically rearranged Doppler output
% files are saved
out_path = [file_path,month,'\ave_output'];

if(~exist(out_path,'dir'))
    mkdir(out_path);
end
% combine Doppler output files from each day
%seq_files_sub = struct([]);
p_f = cell(0);
n_f = cell(0);
e_f = cell(0);
for i_date = 1:length(dates)
    p_sub = cell(0);
    n_sub = cell(0);
    e_sub = cell(0);
    if dates(i_date)<10
        path_to_data=[file_path,month,'\decimated\0',sprintf('%d',dates(i_date)),'\output'];
    else
        path_to_data=[file_path,month,'\decimated\',sprintf('%d',dates(i_date)),'\output'];
    end
    seq_files_sub = dir(fullfile(path_to_data,'*DOPPLER*.mat'));
    for i_sub = 1:length(seq_files_sub)
        [p_sub{i_sub},n_sub{i_sub},e_sub{i_sub}] = fileparts(fullfile(path_to_data,seq_files_sub(i_sub).name));
    end
    p_f = [p_f,p_sub];
    n_f = [n_f,n_sub];
    e_f = [e_f,e_sub];
end

i=1;
j=0;
date = zeros(0);
inc_ang_lin = zeros(0);
azimuth = zeros(0);
corr_ori = zeros(0);
corr_fit = zeros(0);
v_h = zeros(0);
v_c = zeros(0);
s = cell(0);
Q_out = zeros(0);
Q_I = zeros(0);
Q_d = zeros(0);
vr_back = zeros(0);
vz_back = zeros(0);
U_out_w = zeros(0);
U_out_uw = zeros(0);
I_c = zeros(0);
b_I = zeros(0);
b_per = zeros(0);
bw = zeros(0);
b_U = zeros(0);
U_d = zeros(0);
vr_offset = zeros(0);
center = zeros(0);
vp = zeros(0);
%while i<2
while i <=length(n_f)-2;
    % combine the two successive sweeps and calculate the average
    oname1=fullfile(out_path,['ave',n_f{i},e_f{i}]);
    try
        covis1 = load([p_f{i},'\',n_f{i},e_f{i}]);
    catch me
        i=i+1;
        disp(['unable to load ', p_f{i},'\',n_f{i},e_f{i}]);
    end
    try
        covis2 = load([p_f{i+1},'\',n_f{i+1},e_f{i+1}]);
    catch me
        i=i+1;
        disp(['unable to load ', p_f{i+1},'\',n_f{i+1},e_f{i+1}])
        continue
    end
    date1 = unixtime2mat(covis1.covis.sweep.timestamp{1});
    date2 = unixtime2mat(covis2.covis.sweep.timestamp{1});
    if date2 - date1 <1/24
        % find two consecutive sweeps and combine them
        j=j+1;
        i=i+2;
        covis = covis1.covis;
        vr1 = covis1.covis.grid{1}.vr;
        vr2 = covis2.covis.grid{1}.vr;
        I1 = covis1.covis.grid{2}.v;
        I2 = covis2.covis.grid{2}.v;
        v_std1 = covis1.covis.grid{1}.std;
        v_std2 = covis2.covis.grid{1}.std;
        I_std1 = covis1.covis.grid{2}.std;
        I_std2 = covis2.covis.grid{2}.std;
        vr1(isnan(vr1))=0;
        vr2(isnan(vr2))=0;
        I1(isnan(I1))=0;
        I2(isnan(I2))=0;
        v_std1(isnan(v_std1))=0;
        v_std2(isnan(v_std2))=0;
        I_std1(isnan(I_std1))=0;
        I_std2(isnan(I_std2))=0;
        covis.grid{1}.vr = (vr1+vr2)/2;
        covis.grid{1}.std = ((v_std1+v_std2)/2)/sqrt(2);
        covis.grid{2}.v = (I1+I2)/2;
        covis.grid{2}.std = ((I_std1+I_std2)/2)/sqrt(2);
        date(j) = (date1+date2)/2;
        [inc_ang_lin(j),azimuth(:,j),corr_ori(:,:,j),corr_fit(:,:,j),cd_I,v_h(:,:,j),v_c(:,j),vz,vz_std,vr,vr_offset(j),s{j}]=vertical_velocity(covis,z_min,z_max,'velocity');
        v_c_sub = v_c(i_z,j);
        covis.grid{1}.v=vz;
        covis.grid{1}.std=vz_std;
        covis.grid{2}.v_filt = cd_I;
        covis.grid{1}.vr = vr;
        [Q_out(:,:,j),Q_I(:,j),Q_d(:,:,j),U_out_uw(:,:,j), U_d(:,j),I_c(:,j),vr_back(:,j),vz_back(:,j),b_I(:,j),b_per(:,j),b_U(:,j),bw(:,:,j),center(:,:,j),vp(j)]=doppler_volume_vertical(covis,z_min,z_max);
        Q_out(v_c_sub<0,:,j) = nan;
        U_out_w(v_c_sub<0,:,j) = nan;
        U_out_uw(v_c_sub<0,:,j) = nan;
        b_I(v_c_sub<0,j) = nan;
        
        save([oname1],'covis');
    else
        j=j+1;
        i=i+1;
        covis=covis1.covis;
        vr = covis.grid{1}.vr;
        v_std = covis.grid{1}.std;
        I = covis.grid{2}.v;
        I_std = covis.grid{2}.std;
        vr(isnan(vr))=0;
        v_std(isnan(v_std))=0;
        I(isnan(I))=0;
        I_std(isnan(I_std))=0;
        covis.grid{1}.vr = vr;
        covis.grid{1}.std = v_std;
        covis.grid{2}.v = I;
        covis.grid{2}.std = I_std;
        date(j)=date1;
        [inc_ang_lin(j),azimuth(:,j),corr_ori(:,:,j),corr_fit(:,:,j),cd_I,v_h(:,:,j),v_c(:,j),vz,vz_std,vr,vr_offset(j),s{j}]=vertical_velocity(covis,z_min,z_max,'velocity');
        v_c_sub = v_c(i_z,j);
        covis.grid{1}.v=vz;
        covis.grid{1}.std=vz_std;
        covis.grid{2}.v_filt = cd_I;
        covis.grid{1}.vr=vr;
        [Q_out(:,:,j),Q_I(:,j),Q_d(:,:,j),U_out_uw(:,:,j),U_d(:,j),I_c(:,j),vr_back(:,j),vz_back(:,j),b_I(:,j),b_per(:,j),b_U(:,j),bw(:,:,j),center(:,:,j),vp(j)]=doppler_volume_vertical(covis,z_min,z_max);   
        Q_out(v_c_sub<0,:,j) = nan;
        U_out_w(v_c_sub<0,:,j) = nan;
        U_out_uw(v_c_sub<0,:,j) = nan;
        b_I(v_c_sub<0,j) = nan;
        
        
        save([oname1],'covis');
    end
end



% save the output variables
time_series_name = ['time_series_quad_',year,'_',month,'.mat'];
time_series_dir = 'F:\COVIS\covis_doppler_data\time_series';
save(fullfile(time_series_dir,time_series_name),'Q_out','Q_I','Q_d','U_out_uw','U_d','I_c','inc_ang_lin','azimuth','date','s','b_I','b_U','vr_back','vz_back','vr_offset','bw','center','vp');
end