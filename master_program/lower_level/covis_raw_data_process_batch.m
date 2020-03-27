% This is the main program used to process a month worth of raw sonar data
% recorded in Imaging and Diffuse modes of COVIS to generate gridded data
% as Matlab structures

% written by guangyux@uw.edu in Oct 2019



%% set up raw data directories
month = 9; % month in which the raw data was recorded
raw_dir = sprintf('C:/COVIS/Axial/COVIS_data/raw/raw_data_combine/2018/%02d',month); % the directory in this line has to be changed to match the directory under which the raw data is stored
raw = dir(fullfile(raw_dir));
e = 0;
raw_path = cell(0);
for i = 3:length(raw)
    if ~raw(i).isdir
        continue;
    else
        e = e+1;
        raw_path{e} = fullfile(raw_dir,raw(i).name);
    end
end

%% set up gridded data directories
% Diffuse-flow data
grid_dir_diff1 = sprintf('C:/COVIS/Axial/COVIS_data/processed/Diffuse_flow/2018/%02d',month); % the directory in this line has to be changed to where the gridded data will be saved
% Imaging data
grid_dir_imag1 = sprintf('C:/COVIS/Axial/COVIS_data/processed/Imaging/2018/%02d',month); % the directory in this line has to be changed to where the gridded data will be saved
% Bathymetry data
grid_dir_bath1 = sprintf('C:/COVIS/Axial/COVIS_data/processed/Bathy/2018/%02d',month); % the directory in this line has to be changed to where the gridded data will be saved
if ~exist(grid_dir_diff1,'dir')
    mkdir(grid_dir_diff1)
end
if ~exist(grid_dir_imag1,'dir')
    mkdir(grid_dir_imag1)
end
if ~exist(grid_dir_bath1,'dir')
    mkdir(grid_dir_bath1)
end

e = 0;
grid_dir_diff = cell(0);
raw = dir(fullfile(raw_dir));
for i = 3:length(raw)
    if ~raw(i).isdir
        continue;
    else
        e = e+1;
        grid_dir_diff{e} = fullfile(grid_dir_diff1,raw(i).name,'\grid_data');
        if ~exist(grid_dir_diff{e},'dir')
            mkdir(grid_dir_diff{e})
        end
    end
end

e=0;
grid_dir_imag = cell(0);
for i = 3:length(raw)
    if ~raw(i).isdir
        continue;
    else
        e = e+1;
        grid_dir_imag{e} = fullfile(grid_dir_imag1,raw(i).name,'\grid_data');
        if ~exist(grid_dir_imag{e},'dir')
            mkdir(grid_dir_imag{e})
        end
    end
end

e=0;
grid_dir_bath = cell(0);
for i = 3:length(raw)
    if ~raw(i).isdir
        continue;
    else
        e = e+1;
        grid_dir_bath{e} = fullfile(grid_dir_bath1,raw(i).name,'\grid_data');
        if ~exist(grid_dir_bath{e},'dir')
            mkdir(grid_dir_bath{e})
        end
    end
end

%% Main loop
for k = 1:length(raw_path)
    raw_path1 = raw_path{k};
    raw_file = dir(fullfile(raw_path1,'covis*'));
    covis_out = cell(1,length(raw_file));
    save_flag = zeros(1,length(raw_file));
    for i = 1:length(raw_file)
        [~,file1,ext] = fileparts(raw_file(i).name);
        if strcmp(ext,'.gz')
            untar(fullfile(raw_path1,raw_file(i).name),raw_path1);
            swp_name = file1(1:end-4);
        elseif strcmp(ext,'.7z')
            filename = fullfile(raw_path1,raw_file(i).name);
            [status,result] = system(['"C:\Program Files\7-Zip\7z.exe" -y x ' '"' filename '"' ' -o' '"' raw_path1 '"']);
            swp_name = file1;
        else
            save_flag(i)=0;
            continue
        end
        
        % grid data
        swp_type = swp_name(23:25);
        switch swp_type
            case 'ima'
                disp(['now processing imaging file: ',swp_name])
                try
                    covis_out{i} = covis_imaging_sweep(raw_path1,swp_name,0,0);
                    save_flag(i)=1;
                catch
                    disp(['Bad sweep:',swp_name])
                    save_flag(i) = 0;
                    fclose('all');
                    try
                        rmdir(fullfile(raw_path1,swp_name),'s');
                    catch
                        fprintf(sprintf('cannot remove %s \n',swp_name));
                        continue
                    end
                    continue
                end
                fclose('all');
                try
                    rmdir(fullfile(raw_path1,swp_name),'s');
                catch
                    fprintf(sprintf('cannot remove %s \n',swp_name));
                    continue
                end
            case 'dif'
                disp(['now processing diffuse-flow file: ',swp_name])
                try
                    covis_out{i} = covis_diffuse_sweep(raw_path1,swp_name,0,0);
                    save_flag(i)=1;
                catch
                    disp(['Bad sweep:',swp_name])
                    save_flag(i) = 0;
                    fclose('all');
                    try
                        rmdir(fullfile(raw_path1,swp_name),'s');
                    catch
                        fprintf(sprintf('cannot remove %s \n',swp_name));
                        continue
                    end
                    continue
                end
                fclose('all');
                try
                    rmdir(fullfile(raw_path1,swp_name),'s');
                catch
                    fprintf(sprintf('cannot remove %s \n',swp_name));
                    continue
                end
            case 'bat'
                disp(['now processing bathymetry file: ',swp_name])
                try
                    covis_out{i} = covis_bathy_sweep(raw_path1,swp_name,0,0);
                    save_flag(i)=1;
                catch
                    disp(['Bad sweep:',swp_name])
                    save_flag(i) = 0;
                    fclose('all');
                    try
                        rmdir(fullfile(raw_path1,swp_name),'s');
                    catch
                        fprintf(sprintf('cannot remove %s \n',swp_name));
                        continue
                    end
                    continue
                end
                fclose('all');
                try
                    rmdir(fullfile(raw_path1,swp_name),'s');
                catch
                    fprintf(sprintf('cannot remove %s \n',swp_name));
                    continue
                end
            otherwise
                disp(['unrecognized data type:',swp_name])
                save_flag(i) = 0;
                continue
        end
    end
    
    % save gridded data
    for e = 1:length(raw_file)
        if save_flag(e)
            swp_name = covis_out{e}.sweep.name;
            swp_type = swp_name(23:25);
            covis = covis_out{e};
            grid_name = [swp_name,'.mat'];
            switch swp_type
                case 'ima'
                    save(fullfile(grid_dir_imag{k},grid_name),'covis');
                case 'dif'
                    save(fullfile(grid_dir_diff{k},grid_name),'covis');
                case 'bat'
                    save(fullfile(grid_dir_bath{k},grid_name),'covis');
                otherwise
                    continue
            end
        end
    end
    
    diff_file = dir(fullfile(grid_dir_diff{k},'*.mat'));
    if isempty(diff_file)
        rmdir(fullfile(grid_dir_diff{k},'..'),'s');
    end
    imag_file = dir(fullfile(grid_dir_imag{k},'*.mat'));
    if isempty(imag_file)
        rmdir(fullfile(grid_dir_imag{k},'..'),'s');
    end
    bath_file = dir(fullfile(grid_dir_bath{k},'*.mat'));
    if isempty(bath_file)
        rmdir(fullfile(grid_dir_bath{k},'..'),'s');
    end  
end


