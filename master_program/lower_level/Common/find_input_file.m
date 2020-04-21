function [fullpath] = find_input_file( filename )
%
% Returns the full path to a file `filename` located in the Inputs/ directory
% in this package's source tree.
%

if isdeployed

  %
  % Path to Inputs files in the packaged Docker image
  %

  fullpath = fullfile('/home','covis','pycovis-matlab','Inputs',filename);

else

  %
  % Path to Inputs relative to this file
  %

  [folder, name, ext] = fileparts(mfilename('fullpath'));
  fullpath = fullfile(folder,'..','..','Inputs',filename);

end
