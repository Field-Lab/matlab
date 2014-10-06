% generates data structure containing full, cleaned spiking responses to
% spatiotemporal electrical stimulation
%
% 1) extracts and combines spike times from responses stored in elecResp,
% stimulus pulse timing stored in movie file (labview output), and files
% containing spike times for interpulse time regions ("..._IPIspikes.mat")
%
% 2) (if necessary) removes artifactual responses generated "bounce-back"
% axonal stimulation (for example, see 2011-06-24-5)
%
% 3) checks for and removes duplicate detections that can occur in overlapping or
% near-overlapping analysis regions (beginning and end of region examined in elecResp
% analysis, including spikes that occur at about the same time as the
% stimulus artifact and are detected as partial spikes both before (in IPIspikes.mat) and
% after (in elecResp) the artifact)
%
% 4) saves combined, cleaned spike sequences to a file
%
%
% Movie file should be main analysis folder (e.g. .../Analysis/2013-05-28-3/movie013)
%   
% 

clear all

saveToDisk = false;

rawDataPath = '/Data/';

%addpath(genpath('/snle/home/snl-e/matlab-standard/private/Geoff/'))


%%%% parameters for 2011-05-28-3 %%%%
if 0
    BasePath = '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/';
    EstimData = 'data013';
    movieFilePath = BasePath;
    FileName = '013';
    movieChunk = 41; %41 or 44 (checked for data012 and data013)
    dataPiece = '2013-05-28-3';
    
    cellIDs =           [166   138   349   437   784   916]; %should match cells in specified ei file
    recElecs =          [6     13    22    31    47    62];  %main recording electrodes used in elecResp analysis
    Patterns =          [1     2     3     4     5     6];
    
    removeArtSpikes = false;
end

%%%% parameters for 2011-07-14-0 %%%%
% ipiSpikes for neuron 498 may be +/-2 detected spikes (of 500+) (based on difference
% between IPIspikes vs. meanSub_IPIspikes)
if 1
    BasePath = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/';
    EstimData = 'data013';
    movieFilePath = BasePath;
    FileName = '013';
    movieChunk = 30; %30/33 37/40 (have checked IPIspikes vs. meanSub_IPIspikes for m37)
    dataPiece = '2011-07-14-0';
    
    cellIDs =           [48   187  304  407  498  618  783  918];
    recElecs =          [4    14   24   33   37   45   53   61];
    Patterns =          [1    2    3    4    5    6    7    8];

    removeArtSpikes = false;
end


%%%% parameters for 2011-06-24-5 %%%%
if 0
    BasePath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/';
    EstimData = 'data015';
    movieFilePath = BasePath;
    FileName = '015';
    movieChunk = 30; %30 or 37 (14 spikes removed for 30, 276 spikes removed for 37)
    dataPiece = '2011-06-24-5';
    
    cellIDs =           [1    33   183   332   214   784]; %should match cells in specified ei file
    recElecs =          [1    2    13    23    26    50];  %main recording electrodes used in elecResp analysis
    Patterns =          [1    2    3     4     5     6];
    
    removeArtSpikes = true;
    blank_window = [49 68];
    cellWithArtSpikesInd = 3;
    artPattern = 2;
end

%%%%

repLength = 20000; %length of a movie chunk in samples

isiDupRemThresh = 10; %minimum interspike interval (samples) to consider 2 spike detections to be from separate spikes

%%


nCells = length(cellIDs);

%% load the spike times from between elecResps

pathToIPIspikes = [BasePath EstimData '-m' num2str(movieChunk) '-filtered/m' num2str(movieChunk) '_IPIspikes.mat'];


% pathToIPIspikes2 = [BasePath EstimData '-m' num2str(movieChunk) '-filtered/m' num2str(movieChunk) 'meanSub_IPIspikes.mat'];
% 
% test1 = load(pathToIPIspikes);
% test2 = load(pathToIPIspikes2);
% 
% for ii = 1:length(test1.ipiSpikes)
%     testSpikes1 = test1.ipiSpikes(ii).spikeTimes;
%     testSpikes2 = test2.ipiSpikes(ii).spikeTimes;
%     
%     if length(testSpikes1) ~= length(testSpikes2) || any(abs(testSpikes2 - testSpikes1)>1)
%         disp(['different spike times in ipiSpikes vs meanSub ipiSpikes for ' num2str(ii) 'th cell'])
%     else
%         disp('all clear')
%     end
% 
% end

if ~exist(pathToIPIspikes, 'file')
   error(['Unable to locate file ' pathToIPIspikes '.' 10,...
       'Path may be wrong or ipiSpikes file may not yet exist.' 10 10,...
       'To generate a new ipiSpikes file, use ipiSpikeSort_script.m' 10 10,...
       'EXITING'])
else
   load(pathToIPIspikes);
end


%% load timing of pattern applications

NS_GlobalConstants = struct('SamplingFrequency',20000,'ChipAddresses',[30 31],...
    'NumberOfChannelsPerChip',32,'CurrentRanges',[0.066 0.266 1.07 4.25 16.9 67.1 264 1040]);

pattern_application_times = ...
    extract_stim_occurrence_timing(FileName, movieFilePath, length(Patterns), movieChunk);

%% figure out which occurrences fall in this movie chunk

pathToElecResp = [BasePath filesep EstimData filesep];

% assumes standard naming and directory structure for p-o-m files 
occInMovieChunk = cell(1,length(Patterns));
for ii = 1:length(Patterns)
    occInMovieChunk{ii} = [];
    allODirs = dir([pathToElecResp 'p' num2str(Patterns(ii)) '_*']);
    for jj = 1:length(allODirs)
        if ~isempty(dir([pathToElecResp allODirs(jj).name '/*_m' num2str(movieChunk)])) %there's a file in this folder with this movie number
            %get the occurence number from the directory name
            occNum = str2double(allODirs(jj).name(strfind(allODirs(jj).name,'_o')+2:end));
            occInMovieChunk{ii} = [occInMovieChunk{ii} occNum];
        end
    end
    occInMovieChunk{ii} = sort(occInMovieChunk{ii});
end

%% load spike times from elecResp files

% load an elecResp to determine number of repetitions of this movie chunk
% and which amplitude it corresponds to
test = load([pathToElecResp 'elecResp_n' num2str(cellIDs(1)) '_p'...
    num2str(Patterns(1)) '_o' num2str(occInMovieChunk{1}(1)) '.mat']);
nReps = test.elecResp.stimInfo.nPulses(test.elecResp.stimInfo.movieNos==movieChunk);
amplitude = find(test.elecResp.stimInfo.movieNos==movieChunk);

successRates = cell(1,nCells); %store success rates of target cells for each occurrence

% get the spike times from the elecResp files
spikeTimesPrelim = cell(nReps,nCells);
for i = 1:nCells;
    successRates{i} = zeros(length(occInMovieChunk{i}),1);

    for j = 1:length(Patterns)
        occToGet = occInMovieChunk{j};
        
        % get spike times
        spikeTimes = zeros(length(occToGet),nReps);
        for kk = 1:length(occToGet)
            thisER = load([pathToElecResp 'elecResp_n' num2str(cellIDs(i)) '_p'...
                num2str(Patterns(j)) '_o' num2str(occToGet(kk))]);
            if ~thisER.elecResp.analysis.finalized(amplitude)
                disp(['***analysis for movie in ' thisER.elecResp.names.rrs_short_name ' not finalized ***'])
            end
            
            latencies = thisER.elecResp.analysis.latencies{amplitude};
            
            if i == j %pattern targeting this cell
                successRates{i}(kk) = sum(latencies~=0)/length(latencies);
            end
            spikeTimes(kk,~latencies) = NaN;
            spikeTimes(kk,~~latencies) = latencies(~~latencies)-1+pattern_application_times{Patterns(j),2}(kk);
        end
                
        % Store the elecResp spikes in a preliminary cell array.
        spikeTimesPrelim{1,i} = [spikeTimesPrelim{1,i} spikeTimes(~isnan(spikeTimes(:,1)),1)'];
        for k = 2:nReps
            spikeTimesPrelim{k,i} = [spikeTimesPrelim{k,i} spikeTimes(~isnan(spikeTimes(:,k)),k)'];
        end
    end
end

% plot success rates across occurrences for each targetted cell
figure
hold on
lineColors = hsv(nCells);
for ii = 1:nCells
    plot(successRates{ii}, '-', 'color', lineColors(ii,:), 'linewidth', 2)
end
set(gca,'ylim', [0 1])
title([dataPiece filesep EstimData ' m' num2str(movieChunk)])
xlabel('occurrence')
ylabel('response probability')


% %%%%  hack to check inter-pulse spikes without having finished elecResp
% %%%%  analysis
% test = load([pathToElecResp 'elecResp_n' num2str(cellIDs(1)) '_p'...
%     num2str(Patterns(1)) '_o' num2str(occInMovieChunk{1}(1)) '.mat']);
% nReps = test.elecResp.stimInfo.nPulses(test.elecResp.stimInfo.movieNos==movieChunk);
% 
% spikeTimesPrelim = cell(nReps,nCells);
% for ii = 1:nCells
%     for kk = 1:nReps
%         spikeTimesPrelim{kk,ii} = [];
%     end
% end


%% combine spikes times from elecResps with spike times from ipiSpikes

% binary vectors flagging which spike times were extracted from elecResp files
spikesFromElecResp = cell(size(spikeTimesPrelim));

for ii = 1:nCells
    for jj = 1:nReps
        inRepBin = ipiSpikes(ii).spikeTimes >= (jj-1)*repLength & ipiSpikes(ii).spikeTimes < jj*repLength;
        inRepSpikeTimes = ipiSpikes(ii).spikeTimes(inRepBin) - (jj-1)*repLength;
        
        spikesFromElecResp{jj,ii} = true(size(spikeTimesPrelim{jj,ii}));
       
        spikeTimesPrelim{jj,ii} = [spikeTimesPrelim{jj,ii} inRepSpikeTimes];
        spikesFromElecResp{jj,ii} = [spikesFromElecResp{jj,ii} false(size(inRepSpikeTimes))];
        
        [spikeTimesPrelim{jj,ii} sortOrd] = sort(spikeTimesPrelim{jj,ii});
        spikesFromElecResp{jj,ii} = spikesFromElecResp{jj,ii}(sortOrd);
    end
end


%% remove spikes corresponding to artifactual axon stimulation

if removeArtSpikes
    %remove artifactual responses of neuron 183 to pattern 2
    [spikeTimesPrelim(:,cellWithArtSpikesInd) spikes_removed removed_bin] = remove_estim_artifactual_responses(spikeTimesPrelim(:,cellWithArtSpikesInd),...
        pattern_application_times{artPattern,2}, blank_window, 'verbose', true, 'plot_results', true);
    
    for ii = 1:length(spikesFromElecResp(:,cellWithArtSpikesInd)) % loop through reps
        spikesFromElecResp{ii,cellWithArtSpikesInd} = spikesFromElecResp{ii,cellWithArtSpikesInd}(~removed_bin{ii});
    end
    
end

%% remove duplicate spikes (spikes found in both the elecResp analysis and
% interpulse interval analysis)
%
% overview: identifies groups of spike times that occur within
% chosen window of time (isiDupRemThresh).  Then chooses to keep either the
% spike that originated from interpulse spike detection (ipiSpikes) or, if all spikes
% originated from elecResp analysis, arbitrarily keeps the first.


for ii = 1:nCells
    mbSpikes(ii).spikeTimesPrelim = spikeTimesPrelim(:,ii); %#ok<SAGROW> %for later comparison
    mbSpikes(ii).spikeTimes       = spikeTimesPrelim(:,ii); %#ok<SAGROW>
    mbSpikes(ii).fromElecResp     = spikesFromElecResp(:,ii); %#ok<SAGROW>
    
    % Loop through the rows (reps) of spikeTimes...
    for jj = 1:length(mbSpikes(ii).spikeTimes)
        iWhile = 0;

        STs = mbSpikes(ii).spikeTimes{jj};
        ERs = mbSpikes(ii).fromElecResp{jj};
        ISIs = diff(STs);
                
        %there has to be a way to parallelize this...
        while(any(ISIs < isiDupRemThresh))
            iWhile = iWhile+1;
            if iWhile > 10000 %infinite loop???
                error('unexpectedly large number of while loop iterations -- check to make sure it isn''t infinitely looping')
            end
                
            dup1 = find(ISIs < isiDupRemThresh, 1);
            dup = STs - STs(dup1) >= 0 & STs - STs(dup1) < isiDupRemThresh;
            dup(dup1) = false;
            
            toRemove = dup & ERs;
                        
            %remove any subsequent spikes within isiDupRemThresh that originally come
            %from elecResp analysis
            STs(toRemove) = [];
            ERs(toRemove) = [];
            
            % catch any remaining duplicates (should only include spikes caught in
            % elecResp followed by spikes caught in ipiSpikes)
            if ~any(toRemove)
                if ERs(dup1) && ~ERs(find(dup,1))
                    STs(dup1) = [];
                    ERs(dup1) = [];
                elseif ISIs(dup1) == 0 %same spike found twice in interpulse interval
                    STs(dup1) = [];
                    ERs(dup1) = [];
                else
                    keyboard
                end
            end
            
            % get ISIs for remaining spikes
            ISIs = diff(STs);
        end
        
        %find remaining duplicates
        dups = find(ISIs < isiDupRemThresh, 1);
        if ~isempty(dups)
            keyboard
        end
        
        mbSpikes(ii).spikeTimes{jj} = STs; %#ok<SAGROW>
        mbSpikes(ii).fromElecResp{jj} = ERs; %#ok<SAGROW>
    end
    
    if 0 % useful plots for verifying that duplicate removal behaved as expected
        figure
        hold on
        for jj = 1:length(mbSpikes(ii).spikeTimes) %reps
            for kk = 1:length(mbSpikes(ii).spikeTimesPrelim{jj}) %spikes
                plot(mbSpikes(ii).spikeTimesPrelim{jj}(kk), jj*8 + mod(kk,3)*1, 'ko')
                if spikesFromElecResp{jj,ii}(kk)
                    plot(mbSpikes(ii).spikeTimesPrelim{jj}(kk), jj*8 + mod(kk,3)*1, 'r.')
                end
            end
            
            for kk = 1:length(mbSpikes(ii).spikeTimes{jj})
                plot(mbSpikes(ii).spikeTimes{jj}(kk), jj*8 + 4 + mod(kk,3)*1, 'bo')
                if mbSpikes(ii).fromElecResp{jj}(kk)
                    plot(mbSpikes(ii).spikeTimes{jj}(kk), jj*8 + 4 + mod(kk,3)*1, 'm.')
                end
            end
            plot([0 20000], (jj*8-1)*[1 1], 'k-')
        end
    end
end


%% determine which spikes are electrically-elicited from target pulse
% criteria: spike time must be after leading edge of pulse (latency > 2 as defined in
% elecResp file), within 1 ms of pulse, and must originate from an elecResp file (rather than
% ipiSpike)

for ii = 1:length(mbSpikes) % loop through cells
    stimWindows = [pattern_application_times{ii,2}'+1 pattern_application_times{ii,2}'+20];
    
    for jj = 1:length(mbSpikes(ii).spikeTimes) %loop through reps
        STs = mbSpikes(ii).spikeTimes{jj}; %vector of spikes
        elicited = false(size(STs));
        for kk = 1:length(stimWindows) % flag any spikes in this time window that come from an elecResp
            elicited = elicited | (STs > stimWindows(kk,1) & STs <= stimWindows(kk,2) & mbSpikes(ii).fromElecResp{jj});
        end
        mbSpikes(ii).elicited{jj} = elicited; %#ok<SAGROW>
    end
end


%% add spike ids
for ii = 1:length(cellIDs)
    mbSpikes(ii).cellIDVis = cellIDs(ii); %#ok<SAGROW>
end

%% save structure to disk

if saveToDisk
    mbSpikesSavePath = [BasePath EstimData filesep 'm' num2str(movieChunk) '_mbSpikes.mat'];
    save(mbSpikesSavePath, 'mbSpikes');
end


%% plot raster of responses to electrical stimulation, with
% estim-driven spikes in red


nCells = length(cellIDs);
figure
for ii = 1:nCells
        
    axes('position', [0.1 (nCells - ii + 1)/(nCells+1) 0.8 1/(nCells+2)])
    hold on
    for kk = 1:length(mbSpikes(ii).spikeTimes) %loop through reps
        for jj = 1:length(mbSpikes(ii).spikeTimes{kk}) %loop through spikes
            if mbSpikes(ii).elicited{kk}(jj)
                plot(mbSpikes(ii).spikeTimes{kk}(jj)*[1 1], [kk-1 kk], 'r-', 'LineWidth', 1)
            else
                plot(mbSpikes(ii).spikeTimes{kk}(jj)*[1 1], [kk-1 kk], 'k-', 'LineWidth', 1)
            end
        end
    end
    hold off
    
    set(gca, 'yLim', [-5 nReps+5], 'xlim', [0 repLength])
    if ii == nCells
        xlabel('time (s)')
        set(gca, 'ylim', [0 nReps])
    else
        set(gca, 'xtick', [], 'ylim', [0 nReps])
    end
    ylabel(['neuron ' num2str(cellIDs(ii))])
end



%%  collect spike pairs from spikeTimesPrelim with ISIs falling within a particular range (for
%  later investigation of cause) into structure "violationsToCheck"

isiMin = 0; %in samples
isiMax = 20;

vCount = 0;
clear violationsToCheck

for ii = 1:nCells
    spikesToCheck = mbSpikes(ii).spikeTimes;
    % Loop through the rows of spikeTimes...
    for jj = 1:length(spikesToCheck)
        spikeIntervals = diff(spikesToCheck{jj});
        refViols = find(spikeIntervals <= isiMax);
        if any(spikeIntervals(refViols)>=isiMin)
            refViolVals = spikeIntervals(refViols);
            for kk = 1:length(refViolVals)
                if refViolVals(kk)>=isiMin
                    vCount = vCount+1;
                    violationsToCheck(vCount).cellInd = ii;
                    violationsToCheck(vCount).rep = jj;
                    violationsToCheck(vCount).spikeTimes =...
                        spikesToCheck{jj}(find(spikeIntervals==refViolVals(kk))+[0 1]);
                end
            end
        end
    end
end

if vCount
    warndlg('detected some remaining refractory period violations -- see violationsToCheck')
end

%% CHECKING PERFORMANCE: pull up some of the raw data to check accuracy of
% spike times


if 0
    filename_movie = [movieFilePath filesep 'movie' FileName];
    full_path = [rawDataPath filesep dataPiece filesep 'data' FileName];
    totalSamples = getRawDataFileSize(full_path)
    
    
    ChunkData = NS_MovieData(FileName, movieChunk, 'fullPath', filename_movie);
    [~, movieBegin, nRepeats, repeatPeriod, movData] = NS_DecodeMovieDataChunk(ChunkData);
    rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
    
    for kk = 1:length(violationsToCheck)
        repNo = violationsToCheck(kk).rep;
        iCell = violationsToCheck(kk).cellInd;
        
        
        sampStart = repeatPeriod*(repNo-1);
        rawData = int16(rawFile.getData(movieBegin+sampStart, repeatPeriod)');
        rawData(1,:) = [];
        
        % compare neuron 1 spike times to data on electrode 1
        spikeTimes = mbSpikes(iCell).spikeTimes{repNo};
        
        figure(5)
        cla
        set(5, 'position', [50 500 1200 600]); hold on
        
        rawData = rawData - int16(mean(rawData,2)*ones(1,size(rawData,2)));
        
        lineColors = lines(7);
        
        set(gca, 'xlim', [0 200] + 100*floor((violationsToCheck(kk).spikeTimes(1)-50)/100))
        %set(gca, 'xlim', [5000 10000])
        
        clusterElecs = getCluster(recElecs(iCell));
        for ii = 2:length(clusterElecs)
            plot(0:size(rawData,2)-1, rawData(clusterElecs(ii),:), '-', 'color', lineColors(ii,:));
        end
        plot(0:size(rawData,2)-1, rawData(clusterElecs(1),:), '-', 'color', lineColors(1,:), 'linewidth', 2);
        
        for ii = 1:length(spikeTimes)
            plot(spikeTimes(ii), -150, 'r*')
        end
        
        drawnow
        
        ylimHigh = get(gca,'ylim'); ylimHigh = ylimHigh(2);
        for ii = 1:length(Patterns)
            pTimes = pattern_application_times{ii,2};
            for jj = 1:length(pTimes)
                text(pTimes(jj), ylimHigh-10, ['p' num2str(ii) ' o' num2str(occInMovieChunk{ii}(jj))])
            end
        end
        
        uiwait(5)
    end
end


%%





