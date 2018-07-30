
[out_path,out_name] = covis_extract(testfiles('diffuse', 'gz'), test_tempdir());

%% Run sweep
metadata = postproc_metadata();

diffuseMatFile = covis_diffuse_sweep( fullfile(out_path, out_name), test_tempdir(), ...
                                    'json_file', "../../Common/input/covis_diffuse.json", ...
                                    'metadata', metadata);
assert(~isempty(diffuseMatFile), "covis_diffuse_sweep returned an empty .mat file path")

validate_imaging_mat( diffuseMatFile )

imgFile = covis_diffuse_plot(diffuseMatFile, test_tempdir());
assert(~isempty(imgFile), "covis_diffuse_plot returned an empty imgfile")
