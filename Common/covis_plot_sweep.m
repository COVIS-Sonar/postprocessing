function imgfile = covis_plot_sweep(matfile, outputdir, varargin)
%

% Check for other args
% p = inputParser;
% parse(p, varargin{:})

% pick a mat file, if none given
if(isempty(matfile))
  error("Matfile %s not specified")
  return
end

% check that archive dir exists
if(~exist(matfile))
    error('Covis .mat file \"%s\" does not exist', matfile);
    return;
end

% load the covis gridded data
load(matfile);

switch(lower(covis.sweep.mode))

    % imaging mode
    case {'imaging', 'dockimaging'}
        % json_proc_file = fullfile('input','covis_image.json');
        % covis = covis_imaging_sweep_kgb(swp_path, swp_name, json_proc_file);
        imgfile = covis_imaging_plot(covis, outputdir, varargin{:});

        % diffuse mode
    case {'diffuse', 'dockdiffuse', 'sonartest'}
        % json_proc_file = fullfile('input','covis_diffuse.json');
        imgfile = covis_diffuse_plot(covis, outputdir, varargin{:});

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
% Ensure imgfile is a char vector (for now)
imgfile = char(imgfile)

end
