function result = poolJobRemote(sched)
%This submits the entire pooljob to scheduler.  
%Pools can be started entirely on the remote resource, so it's possible to
%just package the code containing the pool and run it remotely, with only
%minimal modification to the source (making sure the number of workers in
%the matlabpool job are the same as the requested workers and ensuring
%output arguments are correctly specified).


%We create a parallel job and set the number of tasks.
%Generally speaking, setting the Max and Min number of workers to the same
%number is probably a good idea.  Since we are running locally, we are
%limited to the number of processors on your machine or 8, which ever is
%smaller
j = createMatlabPoolJob(sched);
set(j, 'MaximumNumberOfWorkers', 4);
set(j, 'MinimumNumberOfWorkers', 4);

%The submitted job will run the pooljob.m function.  This file
%must be accessible to workers, so if this job is being run on TUC, then
%the file would need to moved to our system.  You would then need to ensure
%that the file uploaded to TUC was on your path
%sendFileToCAC('examples\poolJobRemoteWorker.m');
set(j,'PathDependencies',{'\\storage01\matlab\naw47'});
%This creates a single task of job 1 that returns 1 output arguments and
%takes no input arguments.
t=createTask(j, @poolJobRemoteWorker, 1, {}) ;
%Since our function fprints to the command window, we'll capture that
alltasks = get(j, 'Tasks');
set(alltasks, 'CaptureCommandWindowOutput', true);
%Lastly, just submit the job.
submit(j);

%submit is asynchronous, so it returns as soon as the job is submitted, you
%have to check when the job completes.  The simplest way to do that is to
%use the provided waitForState function
waitForState(j,'finished');

a = get(alltasks,'Error')
%Once the job is done, we can look at the results, which should be vectors.
%Again, we must unpack the results ourselves, the result is a cell matrix
%with numberofworkers entries. 
result = getAllOutputArguments(j);

%And display the command windowoutput
a = get(alltasks,'CommandWindowOutput');
disp(a);