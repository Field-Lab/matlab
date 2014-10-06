function datarun = get_frame_times(datarun, monitor_refreshes_per_trigger)
% GET_FRAME_TIMES    Calculate the frame times based on the triggers, detect problems
%
% usage: datarun = get_frame_times(datarun, [monitor_refreshes_per_trigger])
% 
% Determines the frame times, using datarun.interval and calling
% MONITOR_REFRESH_TIMES.  Can accept MONITOR_REFRESHES_PER_TRIGGER or
% defaults to 100.
%
% Also saves a problem_frames vector indicating if there are blocks of 
% frames for which the timing is suspect (caused by a frame refresh 
% failure, which occasionally happens).  If no frames are suspect, simply 
% saves a single FALSE, so good way to use this is: 
%     if any(datarun.stimulus.problem_frames)
%         % Do something to deal with problem frames.  Like cut them out.
%     end
%
% (In general for our experiments there seems to be one block of frames at 
% the end that is not returned here, probably because there is no final 
% trigger after them.)
%
% 2010-02 phli
%

if nargin < 2
    monitor_refreshes_per_trigger = 100;
end


triggers = datarun.triggers;
monitor_refresh_rate = datarun.stimulus.monitor_refresh;
interval = datarun.stimulus.interval;


[monitor_times, problem_times] = monitor_refresh_times(triggers, monitor_refresh_rate, monitor_refreshes_per_trigger);

datarun.stimulus.frame_times = monitor_times(1:interval:end);


problem_frames = problem_times(1:interval:end);
if all(problem_frames(:) == false)
    datarun.stimulus.problem_frames = false;
else
    datarun.stimulus.problem_frames = problem_frames;
end