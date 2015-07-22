% for p=202:449
%     Outputs10Jun(p)=Outputsjun10Freddy(p);
% end
% 
% for p=1:561
%     sumi(p)=nansum(Outputs10Jun(p).tracesInfo.I);
% end

% Two output files, one for same recording and stimulating electrode,
% one for
agreement = [];
totalTrials = []; 
totalHumanNegatives = [];
numFalsePositives = [];
numTrueNegatives = [];
totalHumanPositives = [];
numTruePositives = [];
numFalseNegatives = [];

for p = 1:size(Outputs10Jun,2)
    flag = 1; 
    algorithmOutput = Outputs10Jun(p);
        
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
%     disp(['number of stimAmps: ' num2str(length(elecResp.stimInfo.stimAmps))]); 
    clear humanLat; 
    humanLat = zeros(size(spikes,1), size(spikes,2) + 1); 
    humanLat(:) = NaN;
    if elecResp.stimInfo.stimAmps(1) == elecResp.stimInfo.stimAmps(2)
        incSize = 2;
    else
        incSize = 1;
    end
    for a = 1 :incSize: size(elecResp.analysis.latencies,1)/incSize;
        try
            if elecResp.stimInfo.stimAmps(1) == elecResp.stimInfo.stimAmps(2)
                lats1 = elecResp.analysis.latencies{2*a-1}; 
                lats2 = elecResp.analysis.latencies{2*a}; 
                humanLat(ceil(a/incSize),2:end) = [lats1(2:end); lats2(2:end)]';
            else
                humanLat(a,1:size(elecResp.analysis.latencies{a},1)) = elecResp.analysis.latencies{a};
            end
        catch
            try
                disp('does this ever happen?'); 
                humanLat(a,1:(end-1)) = elecResp.analysis.latencies{a}; 
                humanLat(a,end) = NaN;
            catch
                if elecResp.stimInfo.stimAmps(a) == elecResp.stimInfo.stimAmps(a-1)
                     lats1 = elecResp.analysis.latencies{2*a-1}; 
                     lats2 = elecResp.analysis.latencies{2*a}; 
                end
                keyboard;
                if flag == 1
                    disp(['problem at  p:' num2str(p)]);
                    flag = 0;
                end
            end
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
percentFalsePositive = numFalsePositives/totalTrials;
percentTrueNegative = numTrueNegatives/totalTrials;
percentTruePositive = numTruePositives/totalTrials;
percentFalseNegatives = numFalseNegatives/totalTrials;
all_vals = [percentAgreement; percentTrueNegative; percentFalsePositive;...
     percentTruePositive; percentFalseNegatives];

errors = [1.96*sqrt(1/totalTrials*percentAgreement*(1-percentAgreement));...
    1.96*sqrt(1/totalHumanNegatives*percentTrueNegative*(1-percentTrueNegative));...
    1.96*sqrt(1/totalHumanNegatives*percentFalsePositive*(1-percentFalsePositive));...
    1.96*sqrt(1/totalHumanPositives*percentTruePositive*(1-percentTruePositive));...
    1.96*sqrt(1/totalHumanPositives*percentFalseNegatives*(1-percentFalseNegatives))];

%%

agreement = [];
totalTrials = []; 
totalHumanNegatives = [];
numFalsePositives = [];
numTrueNegatives = [];
totalHumanPositives = [];
numTruePositives = [];
numFalseNegatives = [];

for p = 1:size(OutputsSameElectrode,2)
    flag = 1; 
    algorithmOutput = OutputsSameElectrode(p);
        
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
%     disp(['number of stimAmps: ' num2str(length(elecResp.stimInfo.stimAmps))]); 
    clear humanLat; 
    humanLat = zeros(size(spikes,1), size(spikes,2) + 1); 
    humanLat(:) = NaN;
    if elecResp.stimInfo.stimAmps(1) == elecResp.stimInfo.stimAmps(2)
        incSize = 2;
    else
        incSize = 1;
    end
    for a = 1 :incSize: size(elecResp.analysis.latencies,1)/incSize;
        try
            if elecResp.stimInfo.stimAmps(1) == elecResp.stimInfo.stimAmps(2)
                lats1 = elecResp.analysis.latencies{2*a-1}; 
                lats2 = elecResp.analysis.latencies{2*a}; 
                humanLat(ceil(a/incSize),2:end) = [lats1(2:end); lats2(2:end)]';
            else
                humanLat(a,1:size(elecResp.analysis.latencies{a},1)) = elecResp.analysis.latencies{a};
            end
        catch
            try
                disp('does this ever happen?'); 
                humanLat(a,1:(end-1)) = elecResp.analysis.latencies{a}; 
                humanLat(a,end) = NaN;
            catch
                if elecResp.stimInfo.stimAmps(a) == elecResp.stimInfo.stimAmps(a-1)
                     lats1 = elecResp.analysis.latencies{2*a-1}; 
                     lats2 = elecResp.analysis.latencies{2*a}; 
                end
                keyboard;
                if flag == 1
                    disp(['problem at  p:' num2str(p)]);
                    flag = 0;
                end
            end
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
percentFalsePositive = numFalsePositives/totalTrials;
percentTrueNegative = numTrueNegatives/totalTrials;
percentTruePositive = numTruePositives/totalTrials;
percentFalseNegatives = numFalseNegatives/totalTrials;
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
