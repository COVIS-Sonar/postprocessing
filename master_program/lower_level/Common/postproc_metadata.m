function metadata = postproc_metadata()
%
% Produces a struct containing meta-information about the current COVIS
% software:  version strings, Git tags, etc.
%
% Inputs: n/a
%
% Outputs:
%   metadata: Output struct containing the fields:
%      verstr:  Output from covis_version() function
%      matlab_version:         Output from matlab "version" function
%      postprocessing_gitrev:  Git hash for this repo
%      postprocessing_gittags: Git tags for this repo (if any)


  metadata = struct;
  metadata.verstr = covis_version().version_string;
  metadata.matlab_version = version;

  % We can only query the Git metadata when running from the Git repo in Matlab
  % When compiled into Python code, we need to "bake" static values into
  % the python lib.  See the "tmp/static_git_info.m" rule in "Deploy/Makefile"

  % Default values
  metadata.postprocessing_gitrev='';
  metadata.postprocessing_gittags='';

  if exist('static_git_info')==0
    [status,cmdout] = system('git rev-parse HEAD');
    if status == 0
      metadata.postprocessing_gitrev = cmdout;
    end

    [status,cmdout] = system('git describe --tags');
    if status == 0
      metadata.postprocessing_gittags = cmdout;
    end

  else
    % If not running from within Matlab, expect static values to be included
    gitinfo = static_git_info();

    metadata.postprocessing_gitrev = gitinfo.gitrev;
    metadata.postprocessing_gittags = gitinfo.gittags;
  end

end
