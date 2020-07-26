
%% matlab -nodisplay -nosplash -r "initial_matlab_process()"

addpath(genpath('master_program'))
addpath(genpath('ThirdParty/GSW/'))

%covis = covis_raw_sweep('TestData/','COVIS-20191024T003346-diffuse3.tar.gz',0)

covis = covis_raw_sweep('TestData/2020/04/21/','COVIS-20200421T133002-diffuse1.7z',0)
save('COVIS-20200421T133002-diffuse1.mat', 'covis')

exit()
