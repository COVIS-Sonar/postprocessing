function [hdr, data] = covis_read( file )
% Read Reson R7038 data file.  
% This function returns a structure containing the file header information 
% and data matrix. The data matrix is of size [number_samples,number_chans].
% The data matrix is complex with each data column representing the 
% discrete time series of the baseband recieved signal from each element.
% Each data sample is a quadrature pair (I&Q) with 12 bits of resolution.
% The header structure (hdr) contains all the data read from the file header. 
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
% Version 1.0 - cjones@apl.washington.edu 06/2010
% Version 2.0 - Mike Kenny 08/2018


format = 'ieee-le';
fp = fopen(file,'r',format);
if(fp <= 0) 
    error('Error openning f7038 data file');
end

hdr.serial_number = fread(fp, 1, 'uint64');
hdr.number = fread(fp, 1, 'uint32');
reserved = fread(fp, 1, 'uint16'); 
hdr.total_nelems = fread(fp, 1, 'uint16');
hdr.nsamps = fread(fp, 1, 'uint32');
hdr.nelems = fread(fp, 1, 'uint16');
hdr.first_samp = fread(fp, 1, 'uint32');
hdr.last_samp = fread(fp, 1, 'uint32');
hdr.samp_type = fread(fp, 1, 'uint16');
reserved = fread(fp, 7, 'uint32');

nsamps = hdr.nsamps;
nelems = hdr.nelems;

% read element number list
hdr.element = fread(fp, nelems, 'uint16');

% read data (I&Q pairs) for all elements 
data = fread(fp, [2*nelems, nsamps], 'int16');

fclose(fp);

% If the sample type field is 16, the data value is stored in
% the upper 12 bits rather than the lower 12 bits.
if(hdr.samp_type == 16)
    data = data / 16;
end

if(size(data,2) < nsamps) 
    fprintf('Error reading %d samples\n', nsamps);
end

% make a complex number from real and imag components
data = data(1:2:end-1,:)' + 1i*(data(2:2:end,:)');

end

