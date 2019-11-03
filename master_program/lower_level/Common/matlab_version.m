function verstr = matlab_version()
%
% Returns a string describing the current Matlab version.
% I had unusual problems calling the Matlab "version" function from
% within the python module.
% 
% Useful as a test function
%

    verstr = version;
