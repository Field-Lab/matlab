function prod1 = poolJobRemoteWorker()
%This is a slightly modified version of pooljob to run remotely in a
%matlabpooljob.  This is spawned inside of a matlabpool job which means
%that you end up with number of workers-1 labs doing computation.

%Let's setup a matrix of random numbers that we might do something with
%This creates a matrix of random numbers with 200 rows and 4000 columns
%Notice this data is used by ALL workers and so must be created once and
%shared by them.
somemat = rand(600,10000);

%Open up the pool
%This forces the local pool and starts 2 processes.  Note in the local
%configuration, you are limited to the number of processors (cores) on the
%local machine or 8, whichever is greater.
%Running inside of a matlabpool job configures the pool for you already, so
%you don't do this.  However, the number of labs (pool workers) is actually
%number of workers-1.  Inside of a PoolJob, calling matlabpool open is an error
%as the pool has already been created for you based on the the number of workers.
%matlabpool open 2

%We are now inside the pool, this has no effect on regular commands
%These are run serially:
sums = sum(somemat,1);
one = somemat*somemat';
[r,c] = size(somemat);

%If we want to execute in parallel, we use parfor
%This is automatically parallelized across the workers in the lab.
%Note that spmd is also supported.
tic;
parfor i = 1:r
    prod1(i) = prod(somemat(i,:));
end 
partime=toc;

%The iterator on the loop defines what a unit of "work" is here, and one
%must be careful how that works.  Your first impulse might be to do
%something like this where we operate on submatrices (say 10 rows at a
%time).  You have to be careful about this in parfors, see
%http://www.mathworks.com/access/helpdesk/help/toolbox/distcomp/brdqtjj-1.h
%tml for information on what you can do with slicing variables.
%parfor i=1:10:r-10
%    prod3(i:i+9) = prod(somemat(:,i:i+9),1);
%end

%We can still perform serial for loops as well
tic;
for i=1:r
    prod2(i) = prod(somemat(i,:));
end
sertime = toc;

for i = 1:r
    if prod1(i) ~= prod2(i)
        fprintf('Detected difference between serial and parallel versions\n');
    end
end

%Whenever we finish, we should close the pool.
%matlabpool close

%Data created inside the pool is available outside the loop
fprintf('Parallel time to sum columns is %d seconds\n', partime);
fprintf('Serial time to sum columns is %d seconds\n',sertime);
fprintf('Number of labs is %d\n', numlabs);