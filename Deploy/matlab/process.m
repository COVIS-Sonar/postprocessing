function [output] = process( inputFile, outputFile )

%
% Some portion of the output structure cannot be passed back to Python
% space.  Instead, do some of the processing in Matlab space and return
% a "python-sanitized" data structure.
%

[filepath,name,ext] = fileparts( inputFile );

covis = covis_raw_sweep(filepath,strcat(name,ext),0)
save( outputFile, 'covis')

output = {}
