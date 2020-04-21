

result = postproc_metadata();

disp(result)

assert(isfield(result,'verstr'))
assert(strlength(result.verstr) > 0)

%% The postprocessed metadata should contain the field postprocessing_gitrev
assert(isfield(result,'postprocessing_gitrev'))

%% The postprocessed metadata should contain the field postprocessing_gitrev
assert(isfield(result,'postprocessing_gittags'))

%% postprocessing_gittags should contains "test"
assert(contains(result.postprocessing_gittags,'test'))
