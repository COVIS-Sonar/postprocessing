% This function is used to build the time series of plume properties 
% during the selected time range




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