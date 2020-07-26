function version_data = covis_version()
%
% Returns a string describing the current version of the COVIS post-
% processing code.  Also useful as a test function.
%

    version_data = struct;
    version_data.version_number = "3.2.0";
    version_data.version_string = strcat('COVIS Post-processing v', version_data.version_number);

    version_data;
