% This script is used to compare electrical stimulation results using 
% Gonzalo's spike sorting algorithm to those obtained with human sorting.
% The comparison is done on a trial by trial basis, and reports the cases
% that used the same electrode for stimulating and recording differently
% than cases that used separate electrodes for stimulating and recording.

% Must load two files from Gonzalo's visit in June 2015 that contain the
% sorting output from the alogirhtm. The required variables are called 
% sorted_diffStimRecElecs and sorted_sameStimRecElecs

%% load data
diffRecStim = load(['/Volumes/Lab/Projects/electrical_stim/GM-sorting-'...
    'validation/2015-06-algorithm-results/sorted_diffStimRecElecs']);
sameRecStim = load(['/Volumes/Lab/Projects/electrical_stim/GM-sorting-'...
    'validation/2015-06-algorithm-results/sorted_sameStimRecElecs']); 

sorted_diffStimRecElecs = diffRecStim.sorted_diffStimRecElecs; clear diffRecStim; 
sorted_sameStimRecElecs = sameRecStim.sorted_sameStimRecElecs; clear sameRecStim; 

%% analyze cases with different recording and stimulating electrodes

agreement = [];
totalTrials = []; 
totalHumanNegatives = [];
numFalsePositives = [];
numTrueNegatives = [];
totalHumanPositives = [];
numTruePositives = [];
numFalseNegatives = [];

for p = 1:size(sorted_diffStimRecElecs,2)
    flag = 1; 
    algorithmOutput = sorted_diffStimRecElecs(p);        
    latencies = cell2mat(algorithmOutput.latencies);
    spikes = cell2mat(algorithmOutput.spikes);
    
    % Load elecResp file
    pathname = algorithmOutput.path;
    neuronId = algorithmOutput.neuronInfo.neuronIds; 
    fname = ['elecResp_n' num2str(neuronId) '_p' ...
        num2str(algorithmOutput.stimInfo.patternNo) '.mat'];
    filename = fullfile(pathname,fname);
    temp = load(filename);
    elecResp = temp.elecResp;
    clear humanLat; 
    humanLat = zeros(size(spikes,1), size(spikes,2) + 1); 
    humanLat(:) = NaN;
    if elecResp.stimInfo.stimAmps(1) == elecResp.stimInfo.stimAmps(2)
        incSize = 2;
    else
        incSize = 1;
    end
    for a = 1 :incSize: size(elecResp.analysis.latencies,1)/incSize; 
        if elecResp.stimInfo.stimAmps(1) == elecResp.stimInfo.stimAmps(2)
%             disp('have to collapse conditions that had the same stimulation ampltiude');
            lats1 = elecResp.analysis.latencies{2*a-1};
            lats2 = elecResp.analysis.latencies{2*a};
            humanLat(ceil(a/incSize),2:end) = [lats1(2:end); lats2(2:end)]';
        else
            humanLat(a,1:size(elecResp.analysis.latencies{a},1)) = elecResp.analysis.latencies{a};
        end
    end
    if size(latencies,2) == (size(humanLat,2) - 1)
        humanLat(:,1) = [];
    end
    humanSpikes = humanLat>0;
    
    agreementM =  (humanSpikes == spikes);
    agreement = nansum([agreement sum(agreementM(:))]); 
%     totalTrials = nansum([totalTrials  numel(humanSpikes)]);
     totalTrials = nansum([totalTrials algorithmOutput.tracesInfo.I]);
    totalHumanNegatives = nansum([totalHumanNegatives sum(sum(humanSpikes == 0))]);
    numFalsePositives = nansum([numFalsePositives sum(spikes(humanSpikes == 0))]);
    numTrueNegatives = nansum([numTrueNegatives sum(spikes(humanSpikes == 0) == 0)]);
    totalHumanPositives = nansum([totalHumanPositives sum(sum(humanSpikes == 1))]);
    numTruePositives = nansum([numTruePositives sum(spikes(humanSpikes == 1))]);
    numFalseNegatives = nansum([numFalseNegatives ...
        sum(spikes(humanSpikes == 1) == 0)]);
end
agreement =  sum(agreement(:)); 
percentAgreement = agreement/totalTrials;
percentFalsePositive = numFalsePositives/totalHumanNegatives;
percentTrueNegative = numTrueNegatives/totalHumanNegatives;
percentTruePositive = numTruePositives/totalHumanPositives;
percentFalseNegatives = numFalseNegatives/totalHumanPositives;
all_vals = [percentAgreement; percentTrueNegative; percentFalsePositive;...
     percentTruePositive; percentFalseNegatives];

errors = [1.96*sqrt(1/totalTrials*percentAgreement*(1-percentAgreement));...
    1.96*sqrt(1/totalHumanNegatives*percentTrueNegative*(1-percentTrueNegative));...
    1.96*sqrt(1/totalHumanNegatives*percentFalsePositive*(1-percentFalsePositive));...
    1.96*sqrt(1/totalHumanPositives*percentTruePositive*(1-percentTruePositive));...
    1.96*sqrt(1/totalHumanPositives*percentFalseNegatives*(1-percentFalseNegatives))];

%% analyze cases with the same recording and stimulating electrodes

agreement = [];
totalTrials = []; 
totalHumanNegatives = [];
numFalsePositives = [];
numTrueNegatives = [];
totalHumanPositives = [];
numTruePositives = [];
numFalseNegatives = [];

for p = 1:size(sorted_sameStimRecElecs,2)
     
    algorithmOutput = sorted_sameStimRecElecs(p);
    latencies = cell2mat(algorithmOutput.latencies);
    spikes = cell2mat(algorithmOutput.spikes);
    
    % Load elecResp file
    pathname = algorithmOutput.path;
    neuronId = algorithmOutput.neuronInfo.neuronIds; 
    fname = ['elecResp_n' num2str(neuronId) '_p' ...
        num2str(algorithmOutput.stimInfo.patternNo) '.mat'];
    filename = fullfile(pathname,fname);
    temp = load(filename);
    elecResp = temp.elecResp;
    clear humanLat; 
    humanLat = zeros(size(spikes,1), size(spikes,2) + 1); 
    humanLat(:) = NaN;
    if elecResp.stimInfo.stimAmps(1) == elecResp.stimInfo.stimAmps(2)
        incSize = 2;
    else
        incSize = 1;
    end
    for a = 1 :incSize: size(elecResp.analysis.latencies,1)/incSize;
        if elecResp.stimInfo.stimAmps(1) == elecResp.stimInfo.stimAmps(2)
            lats1 = elecResp.analysis.latencies{2*a-1};
            lats2 = elecResp.analysis.latencies{2*a};
            humanLat(ceil(a/incSize),2:end) = [lats1(2:end); lats2(2:end)]';
        else
            humanLat(a,1:size(elecResp.analysis.latencies{a},1)) = elecResp.analysis.latencies{a};
        end
    end
    
    if size(latencies,2) == (size(humanLat,2) - 1)
        humanLat(:,1) = [];
    end
    humanSpikes = humanLat>0;
    
    agreementM =  (humanSpikes == spikes);
    agreement = nansum([agreement sum(agreementM(:))]); 
%     totalTrials = nansum([totalTrials  numel(humanSpikes)]);
    totalTrials = nansum([totalTrials algorithmOutput.tracesInfo.I]);
    totalHumanNegatives = nansum([totalHumanNegatives sum(sum(humanSpikes == 0))]);
    numFalsePositives = nansum([numFalsePositives sum(spikes(humanSpikes == 0))]);
    numTrueNegatives = nansum([numTrueNegatives sum(spikes(humanSpikes == 0) == 0)]);
    totalHumanPositives = nansum([totalHumanPositives sum(sum(humanSpikes == 1))]);
    numTruePositives = nansum([numTruePositives sum(spikes(humanSpikes == 1))]);
    numFalseNegatives = nansum([numFalseNegatives ...
        sum(spikes(humanSpikes == 1) == 0)]);
end
agreement =  sum(agreement(:)); 
percentAgreement = agreement/totalTrials;
percentFalsePositive = numFalsePositives/totalHumanNegatives;
percentTrueNegative = numTrueNegatives/totalHumanNegatives;
percentTruePositive = numTruePositives/totalHumanPositives;
percentFalseNegatives = numFalseNegatives/totalHumanPositives;
all_valsS = [percentAgreement; percentTrueNegative; percentFalsePositive;...
      percentTruePositive; percentFalseNegatives];
errorsS = [1.96*sqrt(1/totalTrials*percentAgreement*(1-percentAgreement));...
    1.96*sqrt(1/totalHumanNegatives*percentTrueNegative*(1-percentTrueNegative));...
    1.96*sqrt(1/totalHumanNegatives*percentFalsePositive*(1-percentFalsePositive));...
    1.96*sqrt(1/totalHumanPositives*percentTruePositive*(1-percentTruePositive));...
    1.96*sqrt(1/totalHumanPositives*percentFalseNegatives*(1-percentFalseNegatives))];


%% Plots

groupnames = {'% agreement';'% true negative';'% false positive';...
    '% true positive';'% false negative'};
bw_colormap = [100  27.1  0;
    0  23.1  100]/100;
gridstatus = 'y';
bw_legend = {'axonal activation','somatic activation'}; 
figure; 
barweb([all_vals all_valsS], [errors errorsS], [], ...
    groupnames, 'Spike by spike analysis of algorithm v. human', ...
    'error type', 'proportional agreement', bw_colormap, gridstatus, bw_legend);

