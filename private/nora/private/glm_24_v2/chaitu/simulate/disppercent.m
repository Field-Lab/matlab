% disppercent.m
%
%         by: justin gardner
%       date: 10/05/04
%    purpose: display percent done
%             call with -inf for init
%             call with percent to update
%             call with inf to end
%       e.g.:
%
%disppercent(-inf,'Doing stuff');for i =  1:30;mglWaitSecs(0.1);disppercent(i/30);end;elapsedTime = disppercent(inf);
function retval = disppercent(percentdone,mesg)

retval = nan;
% check command line arguments
if ((nargin ~= 1) && (nargin ~= 2))
  help disppercent;
  return
end

% global for disppercent
global gDisppercent;

% if this is an init then remember time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (percentdone == -inf)
  % global to turn off printing
  global gVerbose;
  gDisppercent.verbose = gVerbose;
  if ~gVerbose,return,end
  % set starting time
  gDisppercent.t0 = clock;
  % default to no message
  if (nargin < 2)
    mrDisp(sprintf('00%% (00:00:00)'));
  else
    mrDisp(sprintf('%s 00%% (00:00:00)',mesg));
  end    

% display total time at end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (percentdone == inf)
  if ~gDisppercent.verbose,return,end
  % get elapsed time
  elapsedTime = etime(clock,gDisppercent.t0);
  % separate seconds and milliseconds
  numms = round(1000*(elapsedTime-floor(elapsedTime)));
  numsecs = floor(elapsedTime);
  % if over a minute then display minutes separately
  if numsecs>60
    nummin = floor(numsecs/60);
    numsecs = numsecs-nummin*60;
    % check if over an hour
    if nummin > 60
      numhours = floor(nummin/60);
      nummin = nummin-numhours*60;
      timestr = sprintf('%i hours %i min %i secs %i ms',numhours,nummin,numsecs,numms);
    else
      timestr = sprintf('%i min %i secs %i ms',nummin,numsecs,numms);
    end
  else
    timestr = sprintf('%i secs %i ms',numsecs,numms);
  end
  % display time string
  mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\btook %s\n',timestr));
  retval = elapsedTime;
% otherwise show update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  if ~gDisppercent.verbose,return,end
  % avoid things that will end up dividing by 0
  if (percentdone >= 1)
    percentdone = .99;
  elseif (percentdone <= 0)
    percentdone = 0.01;
  end
  % display percent done and estimated time to end
  if (gDisppercent.percentdone ~= floor(100*percentdone))
    mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b%02i%% (%s)',floor(100*percentdone),disptime(etime(clock,gDisppercent.t0)*(1/percentdone - 1))));
  end
end
% remember current percent done
gDisppercent.percentdone = floor(100*percentdone);

% display time
function retval = disptime(t)


hours = floor(t/(60*60));
minutes = floor((t-hours*60*60)/60);
seconds = floor(t-hours*60*60-minutes*60);

retval = sprintf('%02i:%02i:%02i',hours,minutes,seconds);
