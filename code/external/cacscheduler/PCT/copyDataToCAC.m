function copyDataToCAC(localdata, remoteDataLocation, clusterHost)
%COPYDATATOCAC() - PCT interface function to move data from the local machine
% to TUC.  This is primarily designed for use by PCT Interface
% functions, end users should use the sendFileToCAC, sendDirToCAC, and
% gridFTP functions to move data.
%
% See also sendFileToCAC, sendDirToCAC, gridFTP, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

%fprintf('gridftp(eg) %s to %s:%s\n',localdata,clusterHost,remoteDataLocation);

fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
lFile = java.io.File(localdata);
fm.put(lFile,remoteDataLocation);





