function stimNeuronList = computeAmplitudeModulation(psthStatsFolder, ...
    alternationLogfilePath, outputFolderStatistics, varargin)
%
% Parameters:
%
% Returns:
%

%% Read optional inputs, format the imput folders

global BIN_SIZE
BIN_SIZE = 5;
RESP_THRESHOLD = 0.5;

% Default values
neuronList = [];
startTimeBaseline = 300;
endTimeBaseline = 500;
psthSplitOffset = 50;
endPeak = 150;
pixelsToUm = 7;
gratingsToDiscard = [];
stimNeuronThreshold = 1;

% Checking the optional parameters
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn', ...
        'Unexpected number of arguments');
    throw(err);
end

% Reading the optional input arguments
for kk=1:(nbin/2)
    if ~ischar(varargin{kk*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{kk*2-1})
        case 'neuronlist'
            neuronList = varargin{kk*2};
        case 'pixelstoum'
            pixelsToUm = varargin{kk*2};
        case 'gratingstodiscard'
            gratingsToDiscard = varargin{kk*2};
        case 'psthsplitoffset'
            psthSplitOffset = varargin{kk*2};
        case 'endstimulationpeak'
            endPeak = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end


% Formatting the folders
if psthStatsFolder(end:end)~=filesep
    psthStatsFolder = [psthStatsFolder filesep];
end
if outputFolderStatistics(end:end)~=filesep
    outputFolderStatistics = [outputFolderStatistics filesep];
end

% If output folders don't exist, make them here
if ~exist(outputFolderStatistics,'dir')
    mkdir(outputFolderStatistics);
end

%% Finding all the neuron data

contentsPsthFolder = dir(psthStatsFolder);
neuronNames = struct('name','');
nNeurons = 0;

for kk=1:length(contentsPsthFolder)
    if strfind(contentsPsthFolder(kk).name,'.mat')
        if ~isempty(neuronList)
            cNeuronID = str2double(contentsPsthFolder(kk).name(7:end-4));
            if nnz(neuronList==cNeuronID)>0
                nNeurons = nNeurons + 1;
                neuronNames(nNeurons).name = contentsPsthFolder(kk).name;
                neuronNames(nNeurons).id = cNeuronID;
            end
        else
            nNeurons = nNeurons + 1;
            neuronNames(nNeurons).name = contentsPsthFolder(kk).name;
            neuronNames(nNeurons).id = str2double(contentsPsthFolder(kk).name(7:end-4));
        end
    end
end

%% Loading the gratings order

gratingsMap = readAlternationLogfile(alternationLogfilePath);
gratingsSizes = sort(double(cell2mat(keys(gratingsMap))));
nGratings = length(gratingsSizes);

%% Quantifying amplitude of the response to alternation

alternationData = struct('gratings', [], 'response', [], 'variances', [], ...
    'lin_response', [], 'non_linearity', [], 'transiency', []);
stimNeuronList = [];

for kk=1:nNeurons
    load(fullfile(psthStatsFolder, neuronNames(kk).name));
    
    % One alternationData structure per neuron
    alternationData.gratings = gratingsSizes*pixelsToUm;
    alternationData.response = zeros(size(gratingsSizes));
    alternationData.variances = zeros(size(gratingsSizes));
    alternationData.lin_response = zeros(size(gratingsSizes));
    alternationData.non_linearity = zeros(size(gratingsSizes));
    alternationData.transiency = zeros(size(gratingsSizes));
    
    for ll=1:nGratings
        % Getting experiment IDs and phases corresponding to these IDs 
        % for one grating size
        cExperimentIDs = gratingsMap(gratingsSizes(ll));
        cExperimentPhases = cExperimentIDs(1,:);
        cExperimentIDs = cExperimentIDs(2,:);
        nExperiments = numel(cExperimentIDs);
        
        % We'll store the relevant metrics for the neurons in those 
        % and at the end take the max over all experiments
        allResp = zeros(nExperiments,1);
        allVars = zeros(nExperiments,1);
        allLinComp = zeros(nExperiments,1);
        allNLComp = zeros(nExperiments,1);
        allTrans = zeros(nExperiments,1);
        
        for mm=1:nExperiments
            [cPSTH1, cPSTH2] = splitPSTHInHalf(obj.PSTHs(cExperimentIDs(mm)), psthSplitOffset);

            % Compute peak, background and total number of spikes for each PSTH
            [peak1, backgroundOverPeak1, total1, varPeak1] = analyzeSpikingPatternPSTH(cPSTH1, endPeak);
            [peak2, backgroundOverPeak2, total2, varPeak2] = analyzeSpikingPatternPSTH(cPSTH2, endPeak);

            % Compute firing statistic metrics
            smallTotal = min(total1, total2);
            largeTotal = max(total1, total2);
            R1 = peak1 - backgroundOverPeak2;
            R2 = peak2 - backgroundOverPeak1;
            P1 = peak1 - backgroundOverPeak1;
            P2 = peak2 - backgroundOverPeak2;
            if peak1 > peak2
                transiencyPeak = peak1;
                transiencyBackground = backgroundOverPeak1;
            else
                transiencyPeak = peak2; 
                transiencyBackground = backgroundOverPeak2;
            end
            largePeak = max(P1, P2);
            smallPeak = min(P1, P2);
            
            % Quantify things: response, linear component, non-linear component
            % transiency
            if R1 > R2
                allResp(mm) =R1;
                allVars(mm) = varPeak1/cPSTH1.numberOfPulses;
            else
                allResp(mm) = R2;
                allVars(mm) = varPeak2/cPSTH2.numberOfPulses;
            end
            
            % If no response, then no other quantities
            if allResp(mm)<RESP_THRESHOLD
                allNLComp(mm) = 0;
                allLinComp(mm) = 0;
                allTrans(mm) = 0;
            else
                if largeTotal > 0
                    allLinComp(mm) = (largeTotal - smallTotal)/(largeTotal + smallTotal);
                else
                    allLinComp(mm) = 0;
                end

                if largePeak  > 0
                    allNLComp(mm) = max(0, smallPeak/largePeak);
                else
                    allNLComp(mm) = 0;
                end

                if (transiencyPeak - transiencyBackground) > 0
                    allTrans(mm) = (transiencyPeak - transiencyBackground)/...
                        (transiencyPeak + transiencyBackground);
                else
                    allTrans(mm) = 0;
                end
            end

%             [M1, M2, P1, P2, VM1, VM2] = ...
%                 analyzeSpikingPatternPSTHs(cPSTH1, cPSTH2);
%               
%             if M1 > M2
%                 allResp(mm) =M1;
%                 allVars(mm) = VM1/cPSTH1.numberOfPulses;
%             else
%                 allResp(mm) = M2;
%                 allVars(mm) = VM2/cPSTH2.numberOfPulses;
%             end
%             
%             if allResp(mm)<RESP_THRESHOLD
%                 allNLComp(mm) = 0;
%                 allLinComp(mm) = 0;
%                 allTrans(mm) = 0;
%             else
%                 if M1>M2
%                     allTrans(mm) = 1 - (M1-P1)/(M1+P1);
%                     allLinComp(mm) = (M1-M2)/(M1+M2);
%                 else
%                     allTrans(mm) = 1 - (M2-P2)/(M2+P2);
%                     allLinComp(mm) = (M2-M1)/(M2+M1);
%                 end
%                 if min(P1,P2)>0
%                     allNLComp(mm) = min(P1,P2)/max(P1,P2);
%                 else
%                     allNLComp(mm) = 0;
%                 end
%             end
        end
        
        [alternationData.response(ll), indMax] = max(allResp);
        alternationData.variances(ll) = allVars(indMax);
        alternationData.lin_response(ll) = max(allLinComp);
        alternationData.non_linearity(ll) = max(allNLComp);
        alternationData.transiency(ll) = max(allTrans);
        
    end
    
    % Discarding unwanted gratings
    [~, indToDiscard] = intersect(alternationData.gratings, gratingsToDiscard);
    alternationData.gratings(indToDiscard) = [];
    alternationData.response(indToDiscard) = [];
    alternationData.variances(indToDiscard) = [];
    
    % Checking if neuron was stimulated
    if max(alternationData.response) > stimNeuronThreshold
        stimNeuronList = [stimNeuronList; neuronNames(kk).id];
    end
    
    % Saving the results
    save(fullfile(outputFolderStatistics, [neuronNames(kk).name(1:(end-4)) '_mod.mat']),...
        'alternationData');
end

end % computeAmplitudeModulation

function [PSTH1, PSTH2] = splitPSTHInHalf(PSTH, psthSplitOffset)
% Splits a PSTH object into two halves, with first half starting at time
% psthSplitOffset (specified in ms).

% bin size ASSUMED to be 5 ms
global BIN_SIZE
psth_start_bin = psthSplitOffset/BIN_SIZE;

% Wrap around the PSTH
PSTH.data = [PSTH.data(psth_start_bin+1:end) PSTH.data(1:psth_start_bin)];
PSTH.variances = [PSTH.variances(psth_start_bin+1:end) PSTH.variances(1:psth_start_bin)];
PSTH.subPSTH = [PSTH.subPSTH(:,psth_start_bin+1:end) PSTH.subPSTH(:,1:psth_start_bin)];

% Split
PSTH1.experimentID = PSTH.experimentID;
PSTH1.numberOfPulses = PSTH.numberOfPulses;
PSTH1.data = PSTH.data(:,1:floor(size(PSTH.data,2)/2));
PSTH1.variances = PSTH.data(:,1:floor(size(PSTH.variances,2)/2));
PSTH1.subPSTH = PSTH.data(:,1:floor(size(PSTH.subPSTH,2)/2));

PSTH2.experimentID = PSTH.experimentID;
PSTH2.numberOfPulses = PSTH.numberOfPulses;
PSTH2.data = PSTH.data(:,floor(size(PSTH.data,2)/2+1):end);
PSTH2.variances = PSTH.data(:,floor(size(PSTH.variances,2)/2+1):end);
PSTH2.subPSTH = PSTH.data(:,floor(size(PSTH.subPSTH,2)/2+1):end);

end % splitPSTHInHalf

function [peak, expectedDuringPeakFromBackround, total, varPeak] = analyzeSpikingPatternPSTH(PSTH, endPeak)
% Analyzes the spiking pattern of a PSTH: 
%   - computes the number of spikes during the first 150ms of the PSTH
%   - computes the expected number of spikes over a 150ms period from the
%   background activity
%   - computes the total number of spikes seen over the duration of the
%   PSTH

global BIN_SIZE
if ~exist('endPeak', 'var')
    endPeak = 150; % time of the end of the peak in ms
end

peak = sum(PSTH.data(1:(endPeak/BIN_SIZE)));
expectedDuringPeakFromBackround = sum(PSTH.data((endPeak/BIN_SIZE + 1):end));
% scale expected number of spikes accordingly
expectedDuringPeakFromBackround = expectedDuringPeakFromBackround * ...
    (length(PSTH.data(1:(endPeak/BIN_SIZE)))/length(PSTH.data((endPeak/BIN_SIZE + 1):end)));
total = sum(PSTH.data);
varPeak = sum(PSTH.variances(1:(endPeak/BIN_SIZE)));

end % analyzeSpikingPatternPSTH

% function [M1, M2, P1, P2, VM1, VM2, VP1, VP2] = analyzeSpikingPatternPSTHs(PSTH1, PSTH2)
% 
% pThreshold = 0.01;
% START_TIME_BASELINE = 250;
% END_TIME_BASELINE = 500;
% 
% %% Baseline statistics
% 
% % Getting the baseline statistics - PSTH 1
% [nSpikesBackground1, semBackground1, nTrials1, nBins1] = getNumberOfSpikes(...
%     START_TIME_BASELINE, END_TIME_BASELINE, PSTH1);
% nSpikesExpected1 = nSpikesBackground1/nBins1; % Number of spikes expected per bin
% 
% % Estimate from SEM is biased, hence the correcting factor.
% % Other assumption made is that what happens in a bin is
% % independent from what happens in others, which is an
% % approximation but should be close to true.
% varSpikesExpected1 = nTrials1/(nTrials1-1)*semBackground1*nTrials1/nBins1;
% 
% % Getting the baseline statistics - PSTH 1
% [nSpikesBackground2, semBackground2, nTrials2, nBins2] = getNumberOfSpikes(...
%     START_TIME_BASELINE, END_TIME_BASELINE, PSTH2);
% nSpikesExpected2 = nSpikesBackground2/nBins2; % Number of spikes expected per bin
% 
% % Estimate from SEM is biased, hence the correcting factor.
% % Other assumption made is that what happens in a bin is
% % independent from what happens in others, which is an
% % approximation but should be close to true.
% varSpikesExpected2 = nTrials2/(nTrials2-1)*semBackground2*nTrials2/nBins2;
% 
% %% Finding number of elicited spikes, both histograms
% % Going over all the bins and finding the statistically significant
% % ones, adding their contribution to the number of spikes elicited.
% 
% global BIN_SIZE
% last_bin = START_TIME_BASELINE/BIN_SIZE - 1;
% 
% % M1: PSTH2 peak vs PSTH1 background
% M1 = 0;
% VM1 = 0;
% for mm=1:last_bin
%     nSpikesBin = PSTH2.data(mm);
%     varBin = PSTH2.variances(mm);
% 
%     tValue = abs((nSpikesBin-nSpikesExpected1)/...
%         sqrt(varSpikesExpected1/nTrials1 + varBin/nTrials2));
%     pValue = 1-tcdf(tValue, min(nTrials1,nTrials2)-1);
% 
%     if pValue < pThreshold
%         if nSpikesBin>nSpikesExpected1 % Stimulation
%             M1 = M1 + (nSpikesBin - nSpikesExpected1);
%             VM1 = VM1 + varBin;
%         end
%     end
% end
% 
% % M2: PSTH1 peak vs PSTH2 background
% M2 = 0;
% VM2 = 0;
% for mm=1:last_bin
%     nSpikesBin = PSTH1.data(mm);
%     varBin = PSTH1.variances(mm);
% 
%     tValue = abs((nSpikesBin-nSpikesExpected2)/...
%         sqrt(varSpikesExpected2/nTrials2 + varBin/nTrials1));
%     pValue = 1-tcdf(tValue, min(nTrials1,nTrials2)-1);
% 
%     if pValue < pThreshold
%         if nSpikesBin>nSpikesExpected2 % Stimulation
%             M2 = M2 + (nSpikesBin - nSpikesExpected2);
%             VM2 = VM2 + varBin;
%         end
%     end
% end
% 
% % P1: PSTH2 peak vs PSTH2 background
% P1 = 0;
% VP1 = 0;
% for mm=1:last_bin
%     nSpikesBin = PSTH2.data(mm);
%     varBin = PSTH2.variances(mm);
% 
%     tValue = abs((nSpikesBin-nSpikesExpected2)/...
%         sqrt(varSpikesExpected2/nTrials2 + varBin/nTrials2));
%     pValue = 1-tcdf(tValue, nTrials2-1);
% 
%     if pValue < pThreshold
%         if nSpikesBin>nSpikesExpected2 % Stimulation
%             P1 = P1 + (nSpikesBin - nSpikesExpected2);
%             VP1 = VP1 + varBin;
%         end
%     end
% end
% 
% % P2: PSTH1 peak vs PSTH1 background
% P2 = 0;
% VP2 = 0;
% for mm=1:last_bin
%     nSpikesBin = PSTH1.data(mm);
%     varBin = PSTH1.variances(mm);
% 
%     tValue = abs((nSpikesBin-nSpikesExpected1)/...
%         sqrt(varSpikesExpected1/nTrials1 + varBin/nTrials1));
%     pValue = 1-tcdf(tValue, nTrials1-1);
% 
%     if pValue < pThreshold
%         if nSpikesBin>nSpikesExpected1 % Stimulation
%             P2 = P2 + (nSpikesBin - nSpikesExpected1);
%             VP2 = VP2 + varBin;
%         end
%     end
% end
% 
% end % analyzeSpikingPatternPSTHs
% 
% function [nSpikes, sem, nTrials, nBins] = getNumberOfSpikes(startTime,endTime,PSTH)
% % startTime, endTime specified in ms. 
% % Assuming bin size of 5 ms in the PSTH, Fs = 20,000 Hz
% % Returns the standard error of the mean.
% 
% x = (floor(startTime/5)+1):(floor(endTime/5));
% nSpikes = sum(PSTH.data(x));
% sem = sum(PSTH.variances(x))/PSTH.numberOfPulses; % Square of SEM. Addition because Poisson distr. assumption
% nTrials = PSTH.numberOfPulses;
% nBins = length(x);
% 
% end % getNumberOfSpikes

function [M1, M2, P1, P2, VM1, VM2, VP1, VP2] = analyzeSpikingPatternPSTHs(PSTH1, PSTH2)

[PSTH11, PSTH12] = splitPSTHInHalf(PSTH1, 0);
[PSTH21, PSTH22] = splitPSTHInHalf(PSTH2, 0);

M1 = sum(PSTH21.data) - sum(PSTH12.data);
VM1 = (sum(PSTH21.variances) - sum(PSTH12.variances));
M2 = sum(PSTH11.data) - sum(PSTH22.data);
VM2 = (sum(PSTH11.variances) - sum(PSTH22.variances));
P1 = sum(PSTH21.data) - sum(PSTH22.data);
VP1 = (sum(PSTH21.variances) - sum(PSTH22.variances));
P2 = sum(PSTH11.data) - sum(PSTH12.data);
VP2 = (sum(PSTH11.variances) - sum(PSTH12.variances));

end % analyzeSpikingPatternPSTHs