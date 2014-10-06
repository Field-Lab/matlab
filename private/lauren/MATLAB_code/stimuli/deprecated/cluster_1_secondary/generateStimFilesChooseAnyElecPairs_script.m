
%written in September 2011 to further test the primary-alone artifact

clear all

%cd '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_experiment'

% choose primary electrode and 1 or 2 "secondary" electrodes
% pElec = 6;
% sElecs = [30 40];

pElec = 39;
sElecs = [38 8 62];


relAmps = [0.0001 0.1 0.5 1];

delayInMs = 30;

electrodes = [pElec sElecs];
nElec = length(electrodes);

% make "array": electrode amplitudes, arranged as electrodes x patterns
array = [];

% pairs patterns &&&
for j = 1:length(relAmps)
    for k = 2:nElec
        %both +
        array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        array(k, end) = relAmps(j);
        array(1, end) = 1;

        %primary +, secondary -
        array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        array(k, end) = -1*relAmps(j);
        array(1, end) = 1;
    end
end

% single electrode stimuli &&&

for ii = 1:nElec
    array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
    array(ii, end) = 1;
    
    array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
    array(ii, end) = -1;
end

% generate movie chunks file

delay = delayInMs*20;

nPatterns = size(array, 2);

if nPatterns*delayInMs > 2000
    error('can''t fit all patterns into single movie chunk')
end

nSamples = ceil(nPatterns*2*delay/2000)*2000;

h = msgbox(['Set trigger interval to ' num2str(nSamples/20000) ' seconds']);
boxPos = get(h, 'position');
set(h, 'position', [400 500 boxPos(3) boxPos(4)])


pattern = [randperm(nPatterns) randperm(nPatterns)];
times = 0:delay:delay*(2*nPatterns-1);

Chunk = NS_MovieChunkGenerationForExperiment(times, nSamples, pattern);
MovieChunksFile = [1 Chunk];



keyboard

fid = fopen('chosen_pairs2_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('chosen_pairs2_patterns','wb','ieee-le.l64');
fwrite(fid,array,'double'); %was previously commented out for some reason...
fclose(fid);

fid = fopen('chosen_pairs2_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);


