function out = cac_cancelJobFcn(sched,job)
%CACCANCELJOB() - cancels a job on the TUC cluster at CAC.  The end user
%should not call this function directly, but instead access these functions
%through the cac scheduler object and PCT interface functions.
%
% EXAMPLES
%  cacsched - create sched object
%  j = createJob(sched);
%  t =  createTask(j, @rand, 1, {42});
%  t2 = createTask(j,@rand,1,{25});
%  j.submit(); 
%  cancel(j);
%
% See also cacsched, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

%A job is composed a set of PCT Tasks.  Each task corresponds to a job
%submitted to the the cluster.  So "CancelJob" actually cancels a set of
%"PCT Tasks" == "CAC Jobs".  So we simple extract each task and call
%cancelTask on each one.

%resp = cell(1,job.Tasks.size());

tasks = job.Tasks;
for i = 1:size(tasks)
    cac_cancelTaskFcn(sched,tasks(i));
end


%schedResponse = strcat(resp);
%dummy response appears to be required.
out = true;

