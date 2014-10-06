function [j,a] = cacparsubmit(sched)
%This is a simple example of a parallel job run in matlab

%Create the job from our scheduler and run it on 10 cores
j = createParallelJob(sched);
set(j, 'MaximumNumberOfWorkers', 4);
set(j, 'MinimumNumberOfWorkers', 4);
%upload the .m file that we need to run, this just sums the columns of a
%magic square.  In practice, you wouldn't want to send the file in your
%batch script every time, it's easier to move those files in a seperate
%function.  This simple function just moves a single file up to your home
%directory on TUC, see gridFTP() for more advanced access to our storage
%system.
CACHome = fileparts(which('runTests'));
colSumF = fullfile(CACHome,'examples','cacColSum.m');

userData = sched.UserData;
clusterHost = userData.ClusterHost;
remoteDataLocation = userData.RemoteDataLocation;

n = regexp(remoteDataLocation,'storage01\\matlab\\(?<userid>\w+)\\(?<path>.+)','names');
if ~isempty(n)
	remPath = ['~/' n.path];
    userHome = ['\\storage01\matlab\' n.userid];
else
	remPath = '~';
    userHome = [];
end
remPath = userData.RemotePath;

%sendFileToCAC(colSumF);
copyDataToCAC(colSumF, remPath, clusterHost);
%Create the task and tell it to run cacColSum
%For a parallel job, only 1 task is created and all workers (10 in this
%case) execute the same task.  This particular task takes no input
%arguments and returns none.
t=createTask(j, @cacColSum, 1, {}) ;
%It's often useful to capture the command window output from the task, this
%will capture the output from disp() and fprintf() commands in the task
%function and package them up for you to view.
alltasks = get(j, 'Tasks');
set(alltasks, 'CaptureCommandWindowOutput', true);
%The job needs to know where the task function is.  By default the TUC
%storage is not on the MATLAB path, so we have to explictly add it to the
%path.  The two lines below extract your home storage from the object
%created by cacsched and adds that to the path.
a = get(sched,'ParallelSubmitFcn');
set(j,'PathDependencies',{remoteDataLocation});
%submit the job
submit(j);
%Now we have to wait for the job to finish, waitForState is one way to do
%that.  This call blocks until the job is completed and then retrieves the
%data from the remote scheduler
waitForState(j);
%This call loads the output arguments (1 for each task) into the workspace
%in a variable a.  a{1} corresponds to the output arguments from task 1,
%etc.
a = getAllOutputArguments(j);
destroy(j);