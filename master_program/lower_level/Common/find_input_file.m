function [fullpath] = find_input_file( filename )
%
% Returns the full path to a file `filename` located in the Inputs/ directory
% in this package's source tree.
%

[folder, name, ext] = fileparts(mfilename('fullpath'));
fullpath = fullfile(folder,'..','..','Inputs',filename);
