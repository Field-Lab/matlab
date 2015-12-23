function [j,a] = cacLogSubmit(sched)

j = createJob(sched);
set(j,'FileDependencies',{'cacLog.m'});
t = createTask(j,@cacLog,1,{8});
set(t,'CaptureCommandWindowOutput',true);
submit(j);
%Now we loop doing qpeek and status
type = findtype('distcomp.jobexecutionstate');
states = type.Strings;
numLines = 1;
while 1

    %JobStates not contains an array of string job states.
    %states contains the Matlab list of job states 4 is running, anything
    %beyond is finisherd in someway or another.  

    if strcmp(j.state,'finished')
        %Job has completed, break out of the loop and do something
        break;
    end
    %You could notify yourself of the polling situation if you like.  This
    %can be useful if you're running a large distributed job that you're
    %trying to track.
    %fprintf('%d jobs not yet complete, pausing then rechecking.\n', sum(JS<=4));
    
    %Do a qpeek bef[ore waiting and just print the output.
    %Here we are a little fancy and detect what is new so that we can only
    %display that.  This essentially then does a mirror of the remote
    %logging information on your local display.  The potential for abuse is
    %obvious, have fun!
    s=qpeek(sched,j);
    [r,c] = size(s);
    if r > 0 && r > numLines
        for i=numLines+1:r
            fprintf('%s\n',s(i,:));
        end
        numLines = r;
    end
    %Presumably we do not need to request status every 100 ms, so we just
    %pause for a couple seconds.
    pause(5);
end
waitForState(j);
a = getAllOutputArguments(j);
destroy(j);