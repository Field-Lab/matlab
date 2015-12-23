function [j,a] = cacsubmit(sched)
%This is a simple distributed job submission example.

%This creates a distributed job object which we can attach tasks to.
%Each task will run on a single cpu and not be able to communicate with
%other tasks and each returns a seperate result.
%This the "embarrassingly parallel" case, where the input arguments to each
%task change slightly for each independent task.
j = createJob(sched);
for i = 1:2
    %In a real case, the input arguments would likely change for each task
    %created, unless your doing sort of simulation with a random seed or
    %some such thing.
    createTask(j,@rand,1,{3,3});
end

%Once the job has been configured, we submit it to the scheduler.
submit(j);
%Now we have to wait for the job to finish, waitForState is one way to do
%that.  This call blocks until the job is completed and then retrieves the
%data from the remote scheduler
waitForState(j);
%This call loads the output arguments (1 for each task) into the workspace
%in a variable a.  a{1} corresponds to the output arguments from task 1,
%etc.
a = getAllOutputArguments(j);
destroy(j); % clean-up