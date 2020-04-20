function metadata = postproc_metadata()
%
% Produces a struct containing meta-information about the current COVIS
% software:  version strings, Git tags, etc.
%
%
% Inputs: n/a
%
% Outputs:
%   metadata: Output struct containing the fields:
%      verstr:  Output from covis_version() function
%      gitrev:  Hash for curre version of code
%      gittags: Git tags for current revision of code (if any)


  metadata = struct;
  metadata.verstr = covis_version();
  metadata.matlab_version = version;

  % We can only query the Git metadata when running from the Git repo in Matlab
  % When compiled into Python code, we need to "bake" static values into
  % the python lib.  See the "tmp/static_git_info.m" rule in "Deploy/makefile"

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
    % If not a git repo, use static values
    gitinfo = static_git_info();

    metadata.postprocessing_gitrev = gitinfo.gitrev;
    metadata.postprocessing_gittags = gitinfo.gittags;
  end

end
