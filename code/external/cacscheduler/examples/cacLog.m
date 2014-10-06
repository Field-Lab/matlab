function out = cacLog(num)
%Simple operation that demonstrates the cacschedulerLog and qpeek
%functionality

%Log something write away so we know the job started
cacschedulerLog('Job starting up!');

for i = 1:num
    pause(5*i);
    cacschedulerLog(sprintf('I waited %ds before logging.  yeah!',5*i));
end

cacschedulerLog('Down with task, headed home for meatloaf.');
%Send the same as an output:
out = 'Down with task, headed home for meatloaf.';
