function matfile = covis_process_sweep(filename, outputdir, varargin)

%
% COVIS_PROCESS_SWEEP
%
% covis = process_covis_sweep(swp_path, swp_file)
%
% Function to process a covis sweep data archive.
%
% The input arguments:
%  <swp_path> specifies the path to the data
%  <swp_name> specifies the archive name (.zip or .gz file)
%
% If swp_name is null ([]) or zero, the user is prompted to select
% a zip archive.
% Examples:
%   covis = covis_process_sweep('/data','covis.zip');
%   covos = covis_process_sweep('/data',[]);
%
% The function first unzips the archive into the given path (swp_path),
% if the data directory does not already exist.
% The sweep.json file found in the archive is then parsed to determine
% what kind of sweep to process: 'imaging', 'diffuse', or 'doppler'.
%
% The sweep is then processed accordingly with one of three main functions:
%  covis_imaging_sweep() for an imaging sweep,
%  covis_diffuse_sweep() for a difusse sweep,
%  covis_doppler_sweep())for a doppler sweep. NOTE: DISABLED - this
%       function has been removed as it is not currently reliable on
%       all computers
%
% Processing is performed using the respective JSON parameter files
% (covis_image.json, covis_diffuse.fson, or covis_doppler.json).
%       NOTE: DOPPLER PROCESSING DISABLED
%
% Depending on the type of processing, either a 3D or 2D grid will be
% defined with the sonr data gridded in rectangular coordinates.
%
% The processed data is then plotted using one of three main plot functions:
%  covis_imaging_plot() for an imaging sweep,
%  covis_diffuse_plot() for a difusse sweep,
%  covis_doppler_plot())for a doppler sweep. (DOPPLER PLOTTING DISABLED)
%
% Plotting parameters are defined with the corresponding
% JSON parameter files (covis_image_plot.json, covis_diffuse_plot.fson, or
% covis_doppler_plot.json).
%
% Imaging and doppler sweeps create a 3D grid. Diffuse sweeps create a 2D
% grid with the z value of the grid coords set to a zero (for the bottom).
%
% The results are contined in a single data structure called 'covis',
% which contains the sweep meta data, gidded data, beamformer params,
% ping data, and more.
% The ping data structure (png) is an array with size corresponding
% to the number of pings in the sweep. Each ping structure contains the
% meta data from the .json file corresponding to the ping data file (.bin).
% For example, the pulse length for the tenth ping in the sweep is
% covis.ping(10).hdr.pulse_width.
%
% ----------
% This program is free software distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY. You can redistribute it and/or modify it.
% Any modifications of the original software must be distributed in such a
% manner as to avoid any confusion with the original work.
%
% Please acknowledge the use of this software in any publications arising
% from research that uses it.
%
% ----------
% Version 1.1 - 2/2011,
%    cjones@apl.washington.edu, drj@apl.washington.edu, bemis@rci.rutgers.edu
%
%


close all; % close all open figures and files

% Check for other args
p = inputParser;
addParameter(p,'json_file','',@isstring);
addParameter(p,'metadata',0,@isstruct);
parse(p, varargin{:})

% Extract a COVIS archive, if it hasn't been unpacked already
[swp_dir, swp_name] = covis_extract(filename, '');
swp_path= fullfile(swp_dir, swp_name);

% parse sweep.json file in data archive
json_str = fileread(fullfile(swp_path, 'sweep.json'));
sweep = parse_json(json_str);

% Process and plot the sweep depending on mode
switch(lower(sweep.mode))

    % imaging mode
    case {'imaging', 'dockimaging'}
        % json_proc_file = fullfile('input','covis_image.json');
        % covis = covis_imaging_sweep_kgb(swp_path, swp_name, json_proc_file);
        matfile = covis_imaging_sweep(swp_path, outputdir, varargin{:});

        % diffuse mode
    case {'diffuse', 'dockdiffuse', 'sonartest'}
        % json_proc_file = fullfile('input','covis_diffuse.json');
        matfile = covis_diffuse_sweep(swp_path, outputdir, varargin{:});

        % doppler mode
    case {'doppler', 'dockdoppler'}
        %fprintf('Input file was DOPPLER mode; processing disabled in this version\n')
        % fprintf('Input file was DOPPLER mode\n');
        % % json_proc_file=fullfile('input','covis_doppler.josn');
        % [covis] = covis_doppler_sweep(swp_path,swp_name,json_proc_file);

        % bathy mode
    case {'bathy'}
        % json_proc_file = fullfile('input','covis_bathy.json');
        % covis = covis_bathy_sweep(swp_path, swp_name, json_proc_file);

    otherwise
        disp('Unknown sweep mode.')

end

% Python wrapper doesn't handle strings properly right now.
% Ensure matfile is a char vector (for now)
matfile = char(matfile)

% clean up by deleting the sweep directory
% if(rm_swp_dir)
%     rmdir(swp_dir,'s');
% end

end
