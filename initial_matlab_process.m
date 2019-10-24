
%% Matlab -nodisplay

addpath(genpath('master_program'))
addpath(genpath('ThirdParty/GSW/'))

 covis = covis_raw_sweep('covis-test-data/','COVIS-20191024T003002-diffuse1.tar.gz',0)
 save('output.mat', 'covis')
