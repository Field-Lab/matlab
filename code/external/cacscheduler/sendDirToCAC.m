function sendDirToCAC(lDir,rDir)
% SENDDIRTOCAC recursively copy directory from local system to cac
%
% sendDirToCAC(lDir,rDir)
%
% lDir   local directory
% rDir   remote directory basename
%
% If rDir already exists the contents of lDir are added, overwriting
% existing files as neccessary.
%
% Examples
%
% Send the contents of local directory named 'test' that must exist in the
% current working directory to a directory called 'test' on TUC. The remote
% directory 'test' will be created if necessary:
%
%   sendDirToCAC('test','test');
%
% Send the contents of local directory 'lTest', which must exist in the
% current working directory, to the directory 'rTest' on TUC. The remote
% directory 'rTest' will be created if necessary:
%
%   sendDirToCAC('lTest','rTest');
%
% Add the contents of local directory 'test2' to remote directory 'rTest'
% (created in previous example). A warning that 'rTest' already exists will
% occur, as will warning about other subdirectories in common between
% 'test2' and 'rtest'.  Remote files will be overwritten if they conflict.
%
%   sendDirToCAC('test2','test');
%
% Remote and absolute paths are accepted for the local directory
%
%   sendDirToCAC('c:\projects\matlab\test','test');
%   sendDirToCAC('subdir/subdir','test');
%
% Absolute or relative remote paths can be used but must exist.  See
% gridFTP to create directory structures manually. The following moves the
% contents of local directory '/home/nate/data' to remote directory
% 'test/subdir' ONLY if test already exists:
%
%   sendDirToCAC('/home/nate/data','test/subdir')
%
% See also sendFileToCAC, gridFTP, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing
% Modified by LSF

%% Check arguments
error(nargchk(2,2,nargin,'struct'));
error(nargoutchk(0,0,nargout,'struct'));

%% Get a FileMover
fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();

%% Create the remote directory
fm.createDirectory(rDir);

%% Move files and directory from local to remote directory
% Get list of local files and directories
% For all local files & directories
%   If file
%   then
%      transfer file
%   else if not current or parent
%      transfer directory
%   end

d = dir(lDir);
for i=1:length(d)
  if ~d(i).isdir
    % A file: send it
    moveFile(fm,rDir,fullfile(lDir,d(i).name));
  elseif ~(strcmp(d(i).name,'.') || strcmp(d(i).name,'..'))
    % A sub-directory: transfer it
    sendDirToCAC(...
      fullfile(lDir,d(i).name),...
      fullfile(rDir,d(i).name));
  end
end

function success = moveFile(fm,base,fname)
lFile = java.io.File(fname);
[~,name,ext] = fileparts(fname);
fprintf('Moving file %s to %s\n', fname,fullfile(base,[name ext]));
success = fm.put(lFile,base);
return