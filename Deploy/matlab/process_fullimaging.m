function [output] = process( inputFile, imagingFile )

%
% Some portion of the output structure cannot be passed back to Python
% space.  Instead, do some of the processing in Matlab space and return
% a "python-sanitized" data structure.
%
% The "fullimaging" version returns a two-part struct
% with covis.imaging = {an imaging struct}
% and  covis.diffuse = {a set of diffuse structs}
%

[filepath,name,ext] = fileparts( inputFile );

covis = covis_raw_sweep(filepath,strcat(name,ext),0)

try
  disp(covis)
  imaging = covis.imaging
  save( imagingFile, 'imaging')

  diffuseFile = strrep(imagingFile, "imaging1", "diffuse_fr1")
  diffuse = covis.diffuse
  save( diffuseFile, 'diffuse')

catch ME
    warning(['Exception: ', ME.message, ' at ', ME.stack(1).file, ' line ', int2str(ME.stack(1).line) ])
end

output = {}
