
%% Matlab -nodisplay

addpath(genpath('master_program'))
addpath(genpath('ThirdParty/GSW/'))

%covis = covis_raw_sweep('TestData/','COVIS-20191024T003346-diffuse3.tar.gz',0)

covis = covis_raw_sweep('TestData/','COVIS-20191024T000002-imaging1.tar.gz',0)
save('COVIS-20191024T000002-imaging1.mat', 'covis')
