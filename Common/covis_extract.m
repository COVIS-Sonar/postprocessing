function [swp_path, swp_name] = covis_extract(filename, outputdir)
%
% Extract a COVIS archive file (zip or tar.gz).
%
%
% Inputs:
%   filename - full archive name
%   outputdir - directory to extracte the contents of filename.
%               If left empty (''), extracts same directory as "filename"
%
% Outputs:
%   swp_path - path to extracted directory
%   swp_name - name of extracted directory
%
% If filename is null or zero, the user is prompted to select
% an archive file.
%
% If it has already been extracted in the same location as filename,
% return with swp_path set to the location of the existing directory,
% and swp_name set to the name of the directory.
%
% If the archive is already extracted in outputdir, return with swp_path
% set to outputdir, and swp_name set to the name of the directory.
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
%  Version 1.0 - 10/2011,  cjones@apl.washington.edu
%

swp_path = 0;
swp_name = 0;

% Ensure it's a string before processing...
filename = string(filename);

% pick a sweep archive, if none given
if(isempty(filename))
    error('No filename specified')
    return

%    [swp_name, swp_path] = uigetfile('*.*','Pick a COVIS Sweep Archive');
%    if(swp_name == 0)
%       return;
%    end;
%    filename = fullfile(swp_path, swp_name);
end

if(~exist(filename, 'file'))
    fprintf('COVIS archive %s does not exist\n', filename);
    return;
end

% Check if filename is a directory.
% If so, split it into the final element (the swp_name) and the preceding
% path (the swp_path), to match the output if it _had_ extracted the
% directory.
if(exist(filename, 'dir'))
    % Strip and trailing delimiters if they exist. This will let fileparts
    % correctly separate the swp_name as the last element in the path
    elems = split(filename, filesep);

    swp_name = elems(end);
    swp_path = join( elems(1:end-1), filesep);

    % return with swp_path set to same location as filename
    fprintf('COVIS archive is already extracted, using %s %s\n', swp_path, swp_name);
    return;
end

% get the file name parts
[swp_path, swp_name, ext] = fileparts(filename);

% if no outputdir is given, extract it in the same location as filename
if(isempty(outputdir))
   outputdir = swp_path;
end

% unzip the sweep archive
if(strcmp(ext,'.zip'))

    % unzip file, if it doesn't exist
    if(~exist(fullfile(outputdir, swp_name)))
        fprintf('Extracting archive into %s\n', outputdir);
        fprintf('Unzipping COVIS archive %s ... ', swp_name);
        files = unzip(filename, outputdir);
        fprintf(' done\n');
        % check if zip file has been renamed
        swp_dir = fileparts(files{1});
        if ~strcmp(swp_dir, fullfile(outputdir, swp_name))
            movefile(swp_dir, fullfile(outputdir, swp_name));
        end
    end
    % set new sweep path
    swp_path = outputdir;

elseif(strcmp(ext,'.gz'))

    % untar the sweep archive, if it doen't exist
    [swp_path, swp_name, ext] = fileparts(fullfile(swp_path,swp_name));
    if(~strcmp(ext,'.tar'))
        fprintf('Unknown Covis archive type\n');
        covis = [];
        return;
    end;
    % untar file, if doesn't exist
    if(~exist(fullfile(outputdir,swp_name)))
        fprintf('Extracting archive into %s\n', outputdir);
        fprintf('Untarring COVIS archive %s ... ', swp_name);
        tarfile = char(gunzip(filename, outputdir));
        files = untar(tarfile, outputdir);
        delete(tarfile);
        fprintf(' done\n');
        % check if tar file has been renamed
        swp_dir = fileparts(files{1});
        if ~strcmp(swp_dir, fullfile(outputdir, swp_name))
            movefile(swp_dir, fullfile(outputdir, swp_name));
        end
    end
    % set new sweep path
    swp_path = outputdir;

elseif(strcmp(ext,'.tar'))

    % untar file, if doesn't exist
    if(~exist(fullfile(outputdir,swp_name)))
        fprintf('Untarring COVIS archive %s ... ', filename);
        files = untar(filename, outputdir);
        fprintf(' done\n');
        % check if tar file has been renamed
        swp_dir = fileparts(files{1});
        if ~strcmp(swp_dir, fullfile(outputdir, swp_name))
            movefile(swp_dir, fullfile(outputdir, swp_name));
        end
    end
    % set new sweep path
    swp_path = outputdir;

  elseif(strcmp(ext,'.7z'))

    error("Cannot handle 7z format natively (yet).  Please uncompress manually first.");
    return;


else
    warning("Unknown COVIS archive type: %s . This will probably cause an error.", ext);
    swp_name = [];
    swp_path = [];

    return;
end

% check that archive dir exists
if(~exist(fullfile(swp_path, swp_name)))
    error('Sweep directory does not exist after extraction\n');
    return;
end


end
