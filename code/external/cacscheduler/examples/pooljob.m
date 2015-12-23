function prod2 = pooljob()
%This is a simple example of a matlabpool job someone might have this is
%the basis for a simple tutorial for converting jobs from matlabpools to
%something that will run on TUC.


%Let's setup a matrix of random numbers that we might do something with
%This creates a matrix of random numbers with 200 rows and 4000 columns
%Notice this data is used by ALL workers and so must be created once and
%shared by them.
somemat = rand(600,10000);

%Open up the pool
%This forces the local pool and starts 2 processes.  Note in the local
%configuration, you are limited to the number of processors (cores) on the
%local machine or 8, whichever is greater.
matlabpool open 2
%We are now inside the pool, this has no effect on regular commands
%These are run serially:
sums = sum(somemat,1);
one = somemat*somemat';
[r,c] = size(somemat);

%If we want to execute in parallel, we use parfor
%In this mode, matlab essentially acts as a simple FIFO scheduler and
%throws each piece of work (i, in this case each column of the matrix) at
%the next available processor. 
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

%Note we can still perform serial for loops as well
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
matlabpool close

%Data created inside the pool is available outside the loop
fprintf('Parallel time to sum columns is %d seconds\n', partime);
fprintf('Serial time to sum columns is %d seconds\n',sertime);

