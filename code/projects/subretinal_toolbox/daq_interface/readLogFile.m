function [M, paramNames, allImages, pulseTimes] = readLogFile(logFilePath)
% [M, paramNames, allImages, pulseTimes] = readLogFile(logFilePath)
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
%   - paramNames: a cell array of strings describing what each parameter is
%   - allImages
%
% Note: time is the time of the start of the experiment run, measured in
% seconds elapsed since the operator first hit the run button in the
% Stimulator.
% pulseDuration is the total duration for which the AO was outputting
% samples to the laser or visible led, and returns a value in ms.
%
% Version: v5.01 - 05/29/2013
%


M = [];
paramNames = {'Stimulus','Start Time', 'Experiment duration', 'Stimulus Type', 'Pulse Duration','Number of pulses','Frequency','Power','Image projected'};
allImages = struct('image',[]','h_offset',[],'v_offset',[],'time',[]);
pulseTimes = [];

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
    [M1, param_names1] = parsestimulusstring(stim_str(1:delimiters(1)-2));
    num_trials = cell2mat(M1(:,ismember(param_names1,'Number of trials')));
    freq_stimulus = cell2mat(M1(:,ismember(param_names1,'Frequency')));
    wavelength = cell2mat(M1(:,ismember(param_names1,'Wavelength')));
    start_time = cell2mat(M1(:,ismember(param_names1,'Start time')))*24*3600; % Convert to seconds
    exp_duration = cell2mat(M1(:,ismember(param_names1,'Duration')));

    % Getting the global waveform parameters
    wave_str = stim_str(delimiters(1)+2:delimiters(2)-2);
    [M2, param_names2] = parsewaveformstring(wave_str);
    pulse_times_ms = cell2mat(M2(:,ismember(param_names2,'Pulse times')));
    pulse_width_seq_ms = cell2mat(M2(:,ismember(param_names2,'Pulse duration')));
    irradiance = cell2mat(M2(:,ismember(param_names2,'Frac. of max. irradiance')));  % Test: 2013-03-19
    
    % Cannot put the image data in the matrix M: instead, using a
    % stim_counter, and storing the image information in allImages. 
    pattern_str = stim_str(delimiters(2)+3:end);
    [M3, param_names3] = parsepatternstring(pattern_str);
    image_projected = cell2mat(M3(:,ismember(param_names3,'Images data')));
    
    switch wavelength
        case 'near-IR'
            pulse_type = 0;
        case 'visible'
            pulse_type = 1;
    end
    pulse_duration = pulse_width_seq_ms(end) + pulse_times_ms(end);
    
    
    newRow= [stim_counter start_time exp_duration pulse_type ...
            pulse_duration num_trials freq_stimulus irradiance ...
            stim_counter];
    M = [M; newRow]; %#ok<AGROW>
    allImages(stim_counter) = image_projected;
    pulseTimes = [pulseTimes; pulse_times_ms]; %#ok<AGROW>

    stim_str = fgetl(fid);
    stim_counter = stim_counter + 1;
end
fclose(fid);

% Chaging date reference so that first stimulus starts at time t = 0s
timeExpStart = M(1,ismember(paramNames,'Start Time'));
M(:,ismember(paramNames,'Start Time')) =  M(:,ismember(paramNames,'Start Time'))-timeExpStart;

end % readLogFile