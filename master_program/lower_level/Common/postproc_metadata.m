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

  % We can only query the Git metadata when running from the Git repo
  % when compiled into Python code, we need to "bake" static values into
  % the python lib.  See the "tmp/static_git_info.m" rule in "Deploy/makefile"

  % Default values
  metadata.gitrev='';
  metadata.gittags='';

  if exist('static_git_info')==0
    [status,cmdout] = system('git rev-parse HEAD');
    if status == 0
      metadata.gitrev = cmdout;
    end

    [status,cmdout] = system('git describe --tags');
    if status == 0
      metadata.gittags = cmdout;
    end

  else
    % If not a git repo, try to use static values
    gitinfo = static_git_info();

    metadata.gitrev = gitinfo.gitrev;
    metadata.gittags = gitinfo.gittags;
  end

end
