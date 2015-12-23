function copyDataFromCAC(localdata, remoteDataLocation, clusterHost)
%COPYDATAFROMCAC() - PCT interface function to move data from TUC to the
% local machine.  This is primarily designed for use by PCT Interface
% functions, end users should use the sendFileToCAC, sendDirToCAC, and
% gridFTP functions to move data.
%
% See also sendFileToCAC, sendDirToCAC, gridFTP, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing


%fprintf('gridftp(eg) %s:%s to %s\n',clusterHost,remoteDataLocation,localdata);

%convert localdata to a File object
lFile = java.io.File(localdata);
fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
fm.get(lFile,remoteDataLocation);