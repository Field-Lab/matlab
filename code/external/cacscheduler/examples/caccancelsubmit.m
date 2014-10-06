function j = caccancelsubmit(sched)
%This is a simple distributed job submission example.

%This creates a distributed job object which we can attach tasks to.
%Each task will run on a single cpu and not be able to communicate with
%other tasks and each returns a seperate result.
%for this case we test canceling a task
j = createJob(sched);
for i = 1:2
    % sleep for 10 seconds
    createTask(j,@pause,0,{30});
end

%Once the job has been configured, we submit it to the scheduler.
submit(j);
%Now we have to wait for the job to finish, waitForState is one way to do
%that.  This call blocks until the job is completed and then retrieves the
%data from the remote scheduler
% now we want to pause for a second
pause(2);
% cancel job
cancel(j);
pause(2);
% now wait / cause job to be downloaded
waitForState(j);

