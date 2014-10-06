%function generateStimFilesSpatiotempProbe(centerElecs, timeShiftInMs, periodInMs, offsets)

% centerElecs: vector of 1 or 2 electrodes
% periodInMs: time between application of different stimuli at same center electrode
% timeShiftInMs: offset from time 0 at which to start movie (usually left at 0)
% offsets: vector of time offsets between "prepulse" and stimulating pulse (in ms)
%
% for now, only produces stimuli including 12 nearest-neighboring electrodes for each center, but
% this should be made flexible in the future

cd '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test'

clear all


              
%% parameters that must be set every time                 
           
% first cluster
centerElecs(1) = 60;
amps.centerAmpLims{1} = [0.5 1.4];
centerAmpsThresh(1) = 0.8;

secondaryThreshs{1} = [1     58    54     56; %electrodes
                       1.6   1     1.7    1.8]; %threshold amplitudes
 

% %second cluster
% centerElecs(2) = 53;
% amps.centerAmpLims{2} = [0.7 1.8];
% centerAmpsThresh(2) = 0.1;
% 
% secondaryThreshs{2} = [49     50    51    54; %electrodes
%                        2.1    1.9   2.1   2.3]; %threshold amplitudes

                   
%% parameters that can be left alone

secondaryRange = 0.06; %percent increase and decrease from threshold value for prepulses
neighborAmpsAll = 2; %amplitude of prepulse on all non-primary electrodes in addition to specifically chosen amplitudes
incrementsAll = 0.1; %percent amplitude steps on primary electrode

neighborRadius = 2.1; %in units of nearest-neighbor distance (use 1.8 for 12 electrodes, 2.1 for 18)
periodInMs = 30; %between stimuli on same cluster of electrodes
offsets = [0.5 1 2 4 8]; %in ms

timeShiftInMs = 0; %leave at 0
%% generates specific amplitude values based on chosen parameters
                       
for ii = 1:length(centerElecs)

    amps.neighborAmps{ii} = neighborAmpsAll;
    amps.centerIncrement(ii) = incrementsAll;
    
    amps.centerAmpPrepulses{ii} = [centerAmpsThresh(ii)*(1-secondaryRange) centerAmpsThresh(ii)*(1+secondaryRange) neighborAmpsAll];
    
    specificSecondaryAmps{ii} = [];
    for jj = 1:size(secondaryThreshs{ii},2)
        temp = [secondaryThreshs{ii}(1,jj); (1-secondaryRange)*secondaryThreshs{ii}(2,jj);...
            (1+secondaryRange)*secondaryThreshs{ii}(2,jj); neighborAmpsAll];
        specificSecondaryAmps{ii} = [specificSecondaryAmps{ii}, temp];
    end
end



%%

[electrodes Array details] = generatePatternSpatiotempProbeNew(centerElecs, neighborRadius, amps,...
    'specificSecondaryAmps', specificSecondaryAmps);

keyboard

MovieChunksFile = generateMovieSpatiotempProbeNew(length(centerElecs), timeShiftInMs,...
    periodInMs, offsets, details);


fid = fopen('spatiotemp_probe_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('spatiotemp_probe_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('spatiotemp_probe_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);