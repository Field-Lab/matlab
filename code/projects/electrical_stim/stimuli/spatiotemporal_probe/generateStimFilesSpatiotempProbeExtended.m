
% centerElec: stim pulse electrode (vs prepulse electrode)
% periodInMs: time between application of different stimuli at stim pulse electrode
% offsets: vector of time offsets between "prepulse" and stimulating pulse (in ms)



cd '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_new'

clear all

              
%% parameters that must be set every time                 

% centerElec = 7;
% amps.centerAmpLims = [0.9 3.4];
% centerAmpsThresh = 1.5;
% 
% secondaryThreshs = [6;
%                     3.2];


% centerElec = 39;
% amps.centerAmpLims = [0.7 2.4];
% centerAmpsThresh = 1.1;
% 
% secondaryThreshs = [38    43    42;
%                     2.9   4.4   4.3];

% first cluster
% centerElec = 32;
% amps.centerAmpLims = [0.2 0.45];
% centerAmpsThresh = 0.23;
% 
% secondaryThreshs = [14    22    28    29    33    31; %electrodes
%                     1.0   0.9   0.7   1.2   1.7   1.7]; %threshold amplitudes

% centerElec = 22;
% amps.centerAmpLims = [0.8 1.2];
% centerAmpsThresh = 0.8;
% 
% secondaryThreshs = [];

% centerElec = 39;
% amps.centerAmpLims = [0.4 1];
% centerAmpsThresh = 0.55;
% secondaryThreshs = [35    38    43;
%                     1.9   1.7   1.1];
%  
% centerElec = 24;
% amps.centerAmpLims = [2 5];
% centerAmpsThresh = 2.6;
% secondaryThreshs = [22     21;
%                     2.8    3.5];
%               
% centerElec = 24;
% amps.centerAmpLims = [1.3 3.5];
% centerAmpsThresh = 1.8;
% secondaryThreshs = [22     21;
%                     2      3];
%
% centerElec = 31;
% amps.centerAmpLims = [1.3 4];
% centerAmpsThresh = 2.6;
% secondaryThreshs = [29     30;
%                     3.5     4];

centerElec = 47;
amps.centerAmpLims = [2 6];
centerAmpsThresh = 3.4;
secondaryThreshs = [44      54      31;
                    4.4     4.4     1.8];


                
%% parameters that can be left alone

secondaryRange = 0.06; %percent increase and decrease from threshold value for prepulses
neighborAmpsAll = 4; %amplitude of prepulse on all non-primary electrodes in addition to specifically chosen amplitudes
incrementsAll = 0.2; %percent amplitude steps on primary electrode

%neighborRadius = 2.1; %in units of nearest-neighbor distance (use 1.8 for 12 electrodes, 2.1 for 18)
periodInMs = 128; %between pulses on stim electrode
offsets = [0.5 1 2 4 8 16 32 64]; %in ms

%% generates specific amplitude values based on chosen parameters
                       

amps.neighborAmps = neighborAmpsAll;
amps.centerIncrement = incrementsAll;
    
amps.centerAmpPrepulses = [centerAmpsThresh*(1-secondaryRange) centerAmpsThresh*(1+secondaryRange) neighborAmpsAll];
    
specificSecondaryAmps = [];
for jj = 1:size(secondaryThreshs,2)
    if abs((1+secondaryRange)*secondaryThreshs(2,jj) - neighborAmpsAll)/neighborAmpsAll < 0.05 %higher specified amplitude within 5% of neighborAmpsAll
        temp = [secondaryThreshs(1,jj); (1-secondaryRange)*secondaryThreshs(2,jj); 0; neighborAmpsAll]; %0s won't get used in generateSpatioTempProbeExtended
    else
        temp = [secondaryThreshs(1,jj); (1-secondaryRange)*secondaryThreshs(2,jj);...
            (1+secondaryRange)*secondaryThreshs(2,jj); neighborAmpsAll];
    end
    specificSecondaryAmps = [specificSecondaryAmps, temp]; %#ok<AGROW>
end



%%



[electrodes Array details] = generatePatternSpatiotempProbeExtended(centerElec, amps,...
    'specificSecondaryAmps', specificSecondaryAmps);


MovieChunksFile = generateMovieSpatiotempProbeExtended(periodInMs, offsets, details);

keyboard

fid = fopen('spatiotemp_probe_electrodes2','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('spatiotemp_probe_patterns2','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('spatiotemp_probe_movie2','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);