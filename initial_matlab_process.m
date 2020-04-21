
%% Matlab -nodisplay

addpath(genpath('master_program'))
addpath(genpath('ThirdParty/GSW/'))

%covis = covis_raw_sweep('TestData/','COVIS-20191024T003346-diffuse3.tar.gz',0)

covis = covis_raw_sweep('TestData/2020/04/19/','COVIS-20200419T003001-diffuse1.7z',0)
save('COVIS-20200419T003001-diffuse1.mat', 'covis')
