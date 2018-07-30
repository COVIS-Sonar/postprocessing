function fullpath = input_json_path(filename)

  if getenv('COVIS_IN_DOCKER')
    fullpath = fullfile('input',filename);
  else
    here = mfilename('fullpath');
    [filepath,fname,ext] = fileparts(here);
    fullpath = fullfile(filepath,'input',filename);
  end
