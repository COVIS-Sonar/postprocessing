
%% Matlab -nodisplay

addpath(genpath('master_program'))
addpath(genpath('ThirdParty/GSW/'))

 covis = covis_raw_sweep('TestData/','COVIS-20191024T003346-diffuse3.tar.gz',0)
 save('output.mat', 'covis')
