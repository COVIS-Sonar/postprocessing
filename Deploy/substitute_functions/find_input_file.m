function [fullpath] = find_input_file( filename )
%
% Returns the full path to a file `filename` located in the Inputs/ directory
% in this package's source tree.
%
%
% This version is customized for the Docker image.   It returns an absolute
% path within the image.


fullpath = fullfile('root','pycovis-matlab','Inputs',filename);
