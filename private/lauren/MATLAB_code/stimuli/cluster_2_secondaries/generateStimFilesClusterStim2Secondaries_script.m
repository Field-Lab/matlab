%%%%%%%%% script used to generate stimulus files for spatial stimulation patterns (2 and 3

%%%%%%%%% electrodes)
%%%%%%%%%
%%%%%%%%% code last checked 2010-10-13 by Lauren
%%%%%%%%%
%%%%%%%%% 2011-02-28: added 'smallAmpThresh' to separate small secondary amplitude
%%%%%%%%% patterns from large secondary amplitude patterns

% Variable trigger interval - follow the pop up!

% clear all

%cd '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_experiment'

% choose 1-3 electrodes (can be on edge of the array, but can't be 2nd nearest-neighbors)
%Requirements for clusters of electrodes: must be non-overlapping
centerElectrodes = [47 14]; 
TimeShiftInMs = 0; %offset the first stimulus from the beginning movie chunk
arraySpacing = 60; %30 or 60 micron spacing
includeAnodalPrimary = true; %otherwise, all combinations use primary electrode as cathodal only. setting to true gives 2x as many stimuli
includeTriplets = true; % whether you use patterns with pairs only or pairs + triplets

% For trying extreme ratios
multiplierLimit = 2; %for large relative secondary amplitudes, reduces primary amplitude to keep secondary 
%amplitude within this base amplitude multiplier while maintaining relative secondary:primary amplitude


if includeTriplets
    %%%%%%%%%%%%%%%%%%%%%%%% use this if you want to include patterns
    %%%%%%%%%%%%%%%%%%%%%%%%% involving 2 secondary electrodes
    
    % may want to change vals depending on what you are trying to do. 
    if arraySpacing == 30
        %all relative to the primary
        relAmps.normal = [0.25 0.5 1]; %for 30-micron array % will get applied to all pairs and triplets
        relAmps.pairsOnly = [0.0001 1.5]; %control for primary-alone artifact and test of model limits (very small value used to force rounding down to 0 current)
        relAmps.pairsOnlySmall = 0.1; %control for primary-alone artifact, separated from patterns with larger amplitude pulses on secondary electrodes to avoid rounding errors
    elseif arraySpacing == 60
        relAmps.normal = [0.5 1 2]; %for 60-micron array
        relAmps.pairsOnly = [0.0001 1.5]; %control for primary-alone artifact (very small value used to force rounding down to 0 current)
        relAmps.pairsOnlySmall = 0.1; %control for primary-alone artifact, separated from patterns with larger amplitude pulses on secondary electrodes to avoid rounding errors
    else
        errordlg('invalid array spacing')
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%% use this instead if you only want patterns with
    %%%%%%%%%%%%%%%%%%%%%%%%%%% single secondary electrodes
    if arraySpacing == 30
        relAmps.normal = [];
        relAmps.pairsOnly = [0.0001 0.25 0.5 1 1.5];
        relAmps.pairsOnlySmall = 0.1; %control for primary-alone artifact, separated from patterns with larger amplitude pulses on secondary electrodes to avoid rounding errors
    elseif arraySpacing == 60
        relAmps.normal = [];
        relAmps.pairsOnly = [0.25 0.5 2/3 0.8 1 1.25 1.5 2 4];
        %relAmps.pairsOnlySmall = 0.1; %control for primary-alone artifact, separated from patterns with larger amplitude pulses on secondary electrodes to avoid rounding errors
        relAmps.pairsOnlySmall = [];
    else
        errordlg('invalid array spacing')
    end
end


DelayInMs = 30;


% if any([1 2] == length(centerElectrodes))
%     DelayInMs = 20; %time between repetitions on same primary electrode
% elseif length(centerElectrodes) == 3
%     DelayInMs = 30;
% end

[electrodes Array clusterIDs smallAmpPatterns] = generatePatternClusterStim2SecondariesWrapper(centerElectrodes, relAmps,...
    'includeAnodalPrimaries', includeAnodalPrimary);

for ii = 1:size(Array, 2)
    if any(abs(Array(:,ii))>multiplierLimit);
        Array(:,ii) = Array(:,ii)*(multiplierLimit/max(abs(Array(:,ii))));
    end
end


MovieChunksFile = generateMovieClusterStim2Secondaries(clusterIDs, TimeShiftInMs, DelayInMs, smallAmpPatterns);

keyboard

fid = fopen('2sec_cluster_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('2sec_cluster_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double'); %was previously commented out for some reason...
fclose(fid);

fid = fopen('2sec_cluster_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);


