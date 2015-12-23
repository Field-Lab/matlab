function [times, problem_times] = monitor_refresh_times(triggers, monitor_refresh_rate, monitor_refreshes_per_trigger)
% MONITOR_REFRESH_TIMES    Calculate monitor refresh times based on triggers, detect problems
%
% usage: [times, problem_times] = monitor_refresh_times(triggers, monitor_refresh_rate, monitor_refreshes_per_trigger)
%
% inputs: TRIGGERS                          The trigger times; typically every 100 monitor refreshes
%         MONITOR_REFRESH_RATE              Monitor refresh in Hz; typically in datarun.stimulus.monitor_refresh
%         MONITOR_REFRESHES_PER_TRIGGER     Defaults to 100
%
% outputs: TIMES            The monitor refresh times calculated off triggers
%
%          PROBLEM_TIMES    Same length as TIMES.  Marks with TRUE any
%                           times detected to be problematic (due to
%                           skipped or stuck frames)
%
% 2010-02 phli
%

if nargin < 3
    monitor_refreshes_per_trigger = 100;
end



% Get times of every monitor refresh
times = zeros(monitor_refreshes_per_trigger+1, numel(triggers)-1);
for i = 1:(numel(triggers)-1)
    times(:,i) = linspace(triggers(i), triggers(i+1), monitor_refreshes_per_trigger+1);
end
times = times(1:end-1, :);          % Last row is redundant with first...
times = [times(:); triggers(end)];  % Except for very last point!



if nargout >= 2
    % Check for problematic trigger intervals (dropped/stuck frames)
    trigger_windows = diff(triggers(:));
    monitor_freq = 1 / monitor_refresh_rate;
    calculated_monitor_refreshes_per_trigger = round(trigger_windows ./ monitor_freq);
    problem_trigger_windows = calculated_monitor_refreshes_per_trigger ~= monitor_refreshes_per_trigger;
    
    % If more than half the trigger windows seem wrong, maybe something is off
    if mean(problem_trigger_windows > 0.5)
        warning('MONITOR_REFRESH_TIMES:ManyProblemWindows', 'You have many problem trigger windows; are you sure TRIGGERS and MONITOR_REFRESHES_PER_TRIGGER are correct?');
    end
    
    % Mark all monitor_times that are in problem trigger widows as problems
    % Actually first row should be all false because we can trust the very first element of every trigger window even if frames were dropped within the window...  ToDo?
    problem_times = repmat(problem_trigger_windows', monitor_refreshes_per_trigger, 1);
    problem_times = [problem_times(:); false];
end