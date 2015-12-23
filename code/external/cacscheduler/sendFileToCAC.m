function sendFileToCAC(localFile)
%SENDFILETOCAC(localfile) - send a single local file to TUC storage system
% preserving the name of the file locally, but ignoring any path.  You
% should add a path dependency in your job to refer to this directory:
% set(job,'PathDependencies',{'\\storage01\matlab\username'});
%
% EXAMPLES:
% The following examples all move a single local file (that must exist) to
% your data directory on TUC.
%
% Send a single file in the current directory to TUC
%   sendFileToCAC('colSum.m');
% Send single file in a subdirectory to TUC - note the file is not put in
% subdir/file.mat on TUC!
%   sendFileToCAC('examples/cacMex.m');
% Absolute paths also work:
%   sendFileToCAC('/home/nate/functions/work24.mat');
% 
%  See also sendDirToCAC, gridFTP, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

%Ensure the file exists first
if exist(localFile,'file') == 0
    fprintf('%s doesn''t appear to exist on your filesystem.\nPlease double check the name and path.\n',localFile);
    return
end

fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
lFile = java.io.File(localFile);
success = fm.put(lFile,'~');

if success
    [path,name,ext] = fileparts(localFile); 
    fprintf('%s has been stored on TUC and can be referred to in your \njob at \\\\storage01\\matlab\\USERNAME\\%s\n',localFile,[name ext]);
else
    fprintf('Unable to transfer file.  Check log file for an error.\n');
end