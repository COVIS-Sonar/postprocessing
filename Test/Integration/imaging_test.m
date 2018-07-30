


[out_path,out_name] = covis_extract(testfiles('imaging', 'gz'), test_tempdir());

%% Run sweep
metadata = postproc_metadata();

imagingMatFile = covis_imaging_sweep( fullfile(out_path, out_name), ...
                                    test_tempdir(), ...
                                    'json_file', "../../Common/input/covis_image.json", ...
                                    'metadata', metadata);
assert(~isempty(imagingMatFile), "covis_imaging_sweep returned an empty .mat file path")

validate_imaging_mat( imagingMatFile );

imgFile = covis_imaging_plot(imagingMatFile, test_tempdir());
assert(~isempty(imgFile), "covis_imaging_plot returned an empty imgfile")
