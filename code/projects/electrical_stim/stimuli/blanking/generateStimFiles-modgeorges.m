clear patterns times;

%% 1. General constants - do not modify

numberOfSamples = 10000;
electrodes = 1:512;
array = eye(512, 512);

%% 2. Define of pulse timings and widths

ttlPulseTimings = 10;   % in miliseconds
ttlPulseWidth = 1;      % in miliseconds

laserPulseWidthSequence = [4, 4, 2]; % in milliseconds
laserPulseDelay = [0, 10, 20];       % laser pulse delay after first TTL pulse, in ms
disconnectWidth = 0.8;               % in miliseconds

disconnectTimings = zeros(1, length(laserPulseWidthSequence)*2); % in milliseconds
for kk=1:length(laserPulseWidthSequence);
    disconnectTimings(2*kk - 1) = laserPulseDelay(kk) + ttlPulseTimings(1);
    disconnectTimings(2*kk) = laserPulseDelay(kk) + laserPulseWidthSequence(kk) ...
                                + ttlPulseTimings(1);
end

%% 3. Generate the tables

patterns = [];
times = [];

% 3.1. Define data for TTL pulses
for ii = 1:length(ttlPulseTimings)
    % multiplication by 20 to get value in samples; -1 is the pattern number for TTL
    p1 = ones(1, ttlPulseWidth*20)*(-1); 
    t1 = ttlPulseTimings(ii)*20:(ttlPulseTimings(ii) + ttlPulseWidth)*20 - 1;    
    patterns = [patterns p1];
    times = [times t1];
end

% 3.2. Generate data for disconnection
for ii = 1:length(disconnectTimings)
    % multiplication by 20 to get value in samples; -1 is the pattern number for TTL
    p1 = zeros(1, disconnectWidth*20); 
    t1 = disconnectTimings(ii)*20:(disconnectTimings(ii) + disconnectWidth)*20 - 1;    
    patterns = [patterns p1];
    times = [times t1];
end

%% 4. Generate files - do not modify

chunk = NS_MovieChunkGenerationForExperiment(times, numberOfSamples, patterns);
movieChunksFile = [1 chunk]; % only one movie

fid = fopen(fullfile('.', 'stim_files', 'laser_test1_el'), 'wb', 'l');
fwrite(fid ,electrodes, 'int32');
fclose(fid);

fid = fopen(fullfile('.', 'stim_files', 'laser_test1_pt'), 'wb', 'l');
fwrite(fid, array, 'double');
fclose(fid);

fid = fopen(fullfile('.', 'stim_files', 'laser_test1_mv'), 'wb', 'l');
fwrite(fid, movieChunksFile, 'int32');
fclose(fid);