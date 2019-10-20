% This function is used to process the compressed raw data recorded in 
% various modes of COVIS to generate gridded data and save the output in a 
% Matlab structure array

% version 1.0 by guangyux@uw.edu (Oct 18th, 2019)
%  - based on the code written by Chris Jones in 2010

function covis = covis_raw_sweep(raw_path,raw_name,fig)
% Input:
% raw_path: raw data directory
% raw_name: name of compressed raw data file
% fig: set fig to 1 to plot gridded data. Plotting is muted otherwise.
% Output:
% covis: this is the Matlab structure array that includes the gridded data
% and metadata

% Example
% raw_path = 'raw_path = 'F:\COVIS\Axial\COVIS_data\raw\test';
% raw_name = 'COVIS-20190706T052939-diffuse1.tar.gz';
% fig = 1;

[~,swp_name,ext] = fileparts(raw_name);
if strcmp(ext,'.gz')
    untar(fullfile(raw_path,raw_name),raw_path);
    [~,swp_name,~] = fileparts(swp_name);
elseif strcmp(ext,'.7z')
    filename = fullfile(raw_path,raw_name);
    [~,~] = system(['"C:\Program Files\7-Zip\7z.exe" -y x ' '"' filename '"' ' -o' '"' raw_path '"']);
else
    error('unrecognized data format')
end

% grid data
swp_type = swp_name(23:25);
switch swp_type
    case 'ima'
        disp(['now processing imaging file: ',swp_name])
        try
            covis = covis_imaging_sweep(raw_path,swp_name,0,fig);
        catch
            disp(['Bad sweep:',swp_name])
            try
                rmdir(fullfile(raw_path,swp_name),'s');
            catch
                fprintf(sprintf('cannot remove %s \n',swp_name));
            end
        end
        try
            rmdir(fullfile(raw_path,swp_name),'s');
        catch
            fprintf(sprintf('cannot remove %s \n',swp_name));
        end
    case 'dif'
        disp(['now processing diffuse-flow file: ',swp_name])
        try
            covis = covis_diffuse_sweep(raw_path,swp_name,0,fig);
        catch
            disp(['Bad sweep:',swp_name])
            try
                rmdir(fullfile(raw_path,swp_name),'s');
            catch
                fprintf(sprintf('cannot remove %s \n',swp_name));
            end
        end
        try
            rmdir(fullfile(raw_path,swp_name),'s');
        catch
            fprintf(sprintf('cannot remove %s \n',swp_name));
        end
    case 'bat'
        disp(['now processing bathymetry file: ',swp_name])
        try
            covis = covis_bathy_sweep(raw_path,swp_name,0,fig);
        catch
            disp(['Bad sweep:',swp_name])
            fclose('all');
            try
                rmdir(fullfile(raw_path,swp_name),'s');
            catch
                fprintf(sprintf('cannot remove %s \n',swp_name));
            end
        end
        fclose('all');
        try
            rmdir(fullfile(raw_path,swp_name),'s');
        catch
            fprintf(sprintf('cannot remove %s \n',swp_name));
        end
    otherwise
        disp(['unrecognized data type:',swp_name])
end
end






