function out = cac_cancelTaskFcn(sched,task)
% CACCANCELTASK() - cancels a job on the TUC cluster at CAC.  The end user
% should not call this function directly, but instead access these functions
% through the cac scheduler object and PCT interface functions.
%
% EXAMPLES
%  cacsched - create sched object
%  j = createJob(sched);
%  t =  createTask(j, @rand, 1, {42});
%  t2 = createTask(j,@rand,1,{25});
%  j.submit(); 
%  cancel(t);
%
% See also cacsched, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

%Cancel a "Task" which corresponds to a Scheduler job in this case
%We need to extract the job data from the tasks parent (a job) and then
%index that according to task.id

id = task.ID;
job = task.parent;
jobdata = getJobSchedulerData(sched,job);

if id > size(jobdata.CARObjects,2);
    % more tasks than job IDs
    % skip b/c likely parallel job
    out = true;
    return;
else
    taskdataCarType = jobdata.CARObjects(id);
end

%Now contact the scheduler to cancel the job
hpc = edu.cornell.cac.tuc.matlab.JSDLMediator(1);
schedResponse = char(hpc.jobCancel(taskdataCarType{1}));

%fprintf('HPC Response: %s\n', char(schedResponse));

%dummy response required
out = true;
