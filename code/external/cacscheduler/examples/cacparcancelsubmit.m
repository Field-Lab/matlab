function j = cacparcancelsubmit(sched)
%This is a simple example of a parallel job run in matlab then cancelled

%Create the job from our scheduler and run it on 2 cores
j = createParallelJob(sched);
set(j, 'MaximumNumberOfWorkers', 2);
set(j, 'MinimumNumberOfWorkers', 2);

t=createTask(j, @pause, 0, {30}) ;

%submit the job
submit(j);
%Now we have to wait for the job to finish, waitForState is one way to do
%that.  This call blocks until the job is completed and then retrieves the
%data from the remote scheduler
% now test cancelling
pause(2);
cancel(j);
pause(2);
waitForState(j);
