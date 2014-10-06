function cac_destroyJobFcn(scheduler, job)
%CACDESTROYJOB() - deletes all files related to the provided job locally 
% and on TUC.    The destroy function should be called when all output 
% data has been retrieved from the job. The end user should not call 
% this function directly,  but instead access these functions through 
% the cac scheduler object and PCT interface functions.  
%
% EXAMPLES
%  cacsched - create sched object
%  j = createJob(sched);
%  t =  createTask(j, @rand, 1, {42});
%  t2 = createTask(j,@rand,1,{25});
%  j.submit(); 
%  waitForState(j);
%  a = getAllOutputArguments(j);
%  ONLY AFTER DATA FROM JOB HAS BEEN RETRIEVED
%  destroy(j);
%
% See also cacsched, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing


destroyStart = tic;
userData = getJobSchedulerData(scheduler,job);
clusterHost = scheduler.UserData.ClusterHost;
remoteDataLocation = userData.RemoteDataLocation;
jobLocation = job.pGetEntityLocation;
remPath = userData.RemotePath;

remoteJobDirectory = [ remPath '/' jobLocation ];
remoteJobFiles = [ remPath '/' jobLocation '.*' ];
%fprintf('Destroying Job...\n');
fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
%HACK delete doesn't delete all files in a directory if not empty
%fm.delete([remoteJobDirectory '/matlab_metadata.mat']);
%fm.delete([remoteJobDirectory '/Task*']);
%fm.delete([remoteJobDirectory '/Job*']);
fm.delete([remoteJobDirectory '/*']);
fm.delete(remoteJobDirectory);
fm.delete(remoteJobFiles);
destroyStop = toc(destroyStart);
%fprintf('Destroy took %4.4f seconds\n', destroyStop);
