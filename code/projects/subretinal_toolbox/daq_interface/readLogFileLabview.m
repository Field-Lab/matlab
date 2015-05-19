function M = readLogFileLabview(logFilePath)
% M = readLogFileLabview(logFilePath)
%
% This function reads an experiment text log file and creates a matrix M
% whose columns correspond to the different parameters varied in the
% experiment. 
% Each row of the matrix corresponds to a particular experiment run.
% The function also outputs a cell array with n cells in it, where n is the
% number of columns in n. In each cell there is a string describing what
% the corresponding parameter in M is. 
%
% Input:
%   - logFilePath: a string with the path to the log file formatted to the
%   standard Labview format.
%
% Outputs:
%   - M: a matrix with the experimental parameters used for each run
%
% Note: time is the time of the start of the experiment run, measured in
% seconds elapsed since the operator first hit the run button in the
% Stimulator.
% pulseDuration is the total duration for which the AO was outputting
% samples to the laser or visible led, and returns a value in ms.
%
% Version: v5.01 - 05/29/2013
%


M = struct('stimulus', [], 'start_time', [], 'experiment_duration', [],...
    'stimulus_type', [], 'pulse_durations', [], 'number_of_trials', [],...
    'pulse_times', [], 'frequency', [], 'power', []);

fid = fopen(logFilePath,'r');
if fid<0
    error('Invalid logfile path specified: "%s" cannot be opened.', logFilePath)
end

% Skipping the header
while isempty(strfind(fgetl(fid),'-- Stimuli --'))
end

stim_str = fgetl(fid);
stim_counter = 1;
while stim_str~=-1
    stim_str = stim_str(strfind(stim_str,'::')+3:end);
    delimiters = strfind(stim_str,'::');

    % Getting the global stimulus parameters
    [M1, param_names1] = parsestimulusstringgratings(stim_str(1:delimiters(1)-2));
    num_trials = cell2mat(M1(:,ismember(param_names1,'Number of trials')));
    freq_stimulus = cell2mat(M1(:,ismember(param_names1,'Frequency')));
    wavelength = cell2mat(M1(:,ismember(param_names1,'Wavelength')));
    start_time = cell2mat(M1(:,ismember(param_names1,'Start time')))*24*3600; % Convert to seconds
    exp_duration = cell2mat(M1(:,ismember(param_names1,'Duration')));

    % Getting the global waveform parameters
    wave_str = stim_str(delimiters(1)+2:end);
    [M2, param_names2] = parsewaveformstring(wave_str);
    pulse_times_ms = cell2mat(M2(:,ismember(param_names2,'Pulse times')));
    pulse_width_seq_ms = cell2mat(M2(:,ismember(param_names2,'Pulse duration')));
    irradiance = cell2mat(M2(:,ismember(param_names2,'Frac. of max. irradiance')));  % Test: 2013-03-19
    
    switch wavelength
        case 'near-IR'
            pulse_type = 0;
        case 'visible'
            pulse_type = 1;
    end
    
    % Update 
    M.stimulus(end+1) = stim_counter;
    M.start_time(end+1) = start_time;
    M.experiment_duration(end+1) = exp_duration;
    M.stimulus_type(end+1) = pulse_type;
    M.pulse_durations(end+1,:) = pulse_width_seq_ms;
    M.pulse_times(end+1,:) = pulse_times_ms;
    M.number_of_trials(end+1) =num_trials;
    M.frequency(end+1) = freq_stimulus;
    M.power(end+1,:) = irradiance;

    stim_str = fgetl(fid);
    stim_counter = stim_counter + 1;
end
fclose(fid);

% Changing date reference so that first stimulus starts at time t = 0s
timeExpStart = M.start_time(1);
M.start_time =  M.start_time-timeExpStart;

end % readLogFile