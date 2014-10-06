% script for comparing soma vs. axon thresholds in piece 2008-08-27-2/data002

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/';

noAxonStim = [646 677 858];

midgets = [92 676 766 782 811 886 931];
parasols = [754 902];

midgetAxonPatterns = [18 35 41 23 24 15 6];
parasolAxonPatterns = [27 19];

midgetSomaPatterns = [6 46 53 54 55 60 63];
parasolSomaPatterns = [52 62];

noAxonStimPatterns = [44 46 58];


midgetAxonThresh = zeros(length(midgets), 1);
midgetSomaThresh = zeros(length(midgets), 1);
midgetAxonParams = zeros(length(midgets), 2);
midgetSomaParams = zeros(length(midgets), 2);
midgetAxonMaxResp = zeros(length(midgets), 1);
midgetSomaMaxResp = zeros(length(midgets), 1);

for i = 1:length(midgets)
    temp = load([pathToData filesep 'elecResp_n' num2str(midgets(i)) '_p' num2str(midgetAxonPatterns(i))]);
    elecResp = temp.elecResp;
    %elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', true, 'plotAllFits', true);
    %save([pathToData filesep 'elecResp_n' num2str(midgets(i)) '_p' num2str(midgetAxonPatterns(i))], 'elecResp')
    midgetAxonThresh(i) = elecResp.analysis.threshold;
    midgetAxonParams(i,:) = elecResp.analysis.erfParams;
    midgetAxonMaxResp(i) = max(elecResp.analysis.successRates);
    clear elecResp
    
    temp = load([pathToData filesep 'elecResp_n' num2str(midgets(i)) '_p' num2str(midgetSomaPatterns(i))]);
    elecResp = temp.elecResp;
    %elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', true, 'plotAllFits', true);
    %save([pathToData filesep 'elecResp_n' num2str(midgets(i)) '_p' num2str(midgetSomaPatterns(i))], 'elecResp')
    midgetSomaThresh(i,:) = elecResp.analysis.threshold;
    midgetSomaParams(i,:) = elecResp.analysis.erfParams;
    midgetSomaMaxResp(i) = max(elecResp.analysis.successRates);
    clear elecResp
end

parasolAxonThresh = zeros(length(parasols), 1);
parasolSomaThresh = zeros(length(parasols), 1);
parasolAxonParams = zeros(length(parasols), 2);
parasolSomaParams = zeros(length(parasols), 2);
parasolAxonMaxResp = zeros(length(parasols), 1);
parasolSomaMaxResp = zeros(length(parasols), 1);

for i = 1:length(parasols)
    temp = load([pathToData filesep 'elecResp_n' num2str(parasols(i)) '_p' num2str(parasolAxonPatterns(i))]);
    elecResp = temp.elecResp;
    %elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', true, 'plotAllFits', true);
    %save([pathToData filesep 'elecResp_n' num2str(parasols(i)) '_p' num2str(parasolAxonPatterns(i))], 'elecResp')
    parasolAxonThresh(i) = elecResp.analysis.threshold;
    parasolAxonParams(i,:) = elecResp.analysis.erfParams;
    parasolAxonMaxResp(i) = max(elecResp.analysis.successRates);
    clear elecResp
    
    temp = load([pathToData filesep 'elecResp_n' num2str(parasols(i)) '_p' num2str(parasolSomaPatterns(i))]);
    elecResp = temp.elecResp;
    %elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', true, 'plotAllFits', true);
    %save([pathToData filesep 'elecResp_n' num2str(parasols(i)) '_p' num2str(parasolSomaPatterns(i))], 'elecResp')
    parasolSomaThresh(i) = elecResp.analysis.threshold;
    parasolSomaParams(i,:) = elecResp.analysis.erfParams;
    parasolSomaMaxResp(i) = max(elecResp.analysis.successRates);
    clear elecResp
end

noAxonStimThresh = zeros(length(noAxonStim), 1);
noAxonStimParams = zeros(length(noAxonStim), 2);
noAxonStimMaxResp = zeros(length(noAxonStim), 1);
for i = 1:length(noAxonStim)
    temp = load([pathToData filesep 'elecResp_n' num2str(noAxonStim(i)) '_p' num2str(noAxonStimPatterns(i))]);
    elecResp = temp.elecResp;
    %elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', true, 'plotAllFits', true);
    %save([pathToData filesep 'elecResp_n' num2str(noAxonStim(i)) '_p' num2str(noAxonStimPatterns(i))], 'elecResp')
    noAxonStimThresh(i) = elecResp.analysis.threshold;
    noAxonStimParams(i,:) = elecResp.analysis.erfParams;
    noAxonStimMaxResp(i) = max(elecResp.analysis.successRates);
    clear elecResp
end

allSomasThresh = [midgetSomaThresh; parasolSomaThresh; noAxonStimThresh];
allAxonsThresh = [midgetAxonThresh; parasolAxonThresh; zeros(length(noAxonStimThresh), 1)];
allSomasParams = [midgetSomaParams; parasolSomaParams; noAxonStimParams];
allAxonsParams = [midgetAxonParams; parasolAxonParams];

allSomasMaxResp = [midgetSomaMaxResp; parasolSomaMaxResp; noAxonStimMaxResp];
allAxonsMaxResp = [midgetAxonMaxResp; parasolAxonMaxResp];


stimmdAxonsThresh = [midgetAxonThresh; parasolAxonThresh];
[x somaThreshOrder] = sort(allSomasThresh);
%% plotting: connected thresholds for each cell
%****should have response porbability criterion for axons too, but I know
%that there aren't any axon curves that don't get to >0.5 response
%probability in this dataset

respProbThresh = 0.4; %threshold response probability for cell's data to be plotted

figure
axes('units', 'pixels', 'position', [20 20 161.8*2 200])
hold on
includedSomasThresh = [];
includedAxonsThresh = [];
for i = 1:length(allSomasThresh)
    if allSomasMaxResp(i) >= respProbThresh
        plot(allSomasThresh(somaThreshOrder(i)), i, 'ko', 'markerFaceColor', 'k')
        includedSomasThresh = [includedSomasThresh allSomasThresh(somaThreshOrder(i))];
        if allAxonsThresh(somaThreshOrder(i)) ~= 0
            plot([allSomasThresh(somaThreshOrder(i)) allAxonsThresh(somaThreshOrder(i))], [i i], 'k-')
            plot(allAxonsThresh(somaThreshOrder(i)), i, 'ro',  'markerFaceColor', 'r')
            includedAxonsThresh = [includedAxonsThresh allAxonsThresh(somaThreshOrder(i))];
        end
    else
        disp(['excluded a soma because it didn''t reach ' num2str(respProbThresh) ' probability within analyzed region'])
    end
end
   
% plot([mean(allSomasThresh) + std(allSomasThresh), mean(allSomasThresh) - std(allSomasThresh)], [-3 -3], 'k', 'LineWidth', 2)
% plot(mean(allSomasThresh), -3, 'ks', 'markerFaceColor', 'k')
% plot([mean(stimmdAxonsThresh) + std(stimmdAxonsThresh), mean(stimmdAxonsThresh) - std(stimmdAxonsThresh)], [-2 -2], 'r', 'LineWidth', 2)
% plot(mean(stimmdAxonsThresh), -2, 'rs', 'markerFaceColor', 'r')

plot([mean(includedSomasThresh) + std(includedSomasThresh), mean(includedSomasThresh) - std(includedSomasThresh)], [-3 -3], 'k', 'LineWidth', 2)
plot(mean(includedSomasThresh), -3, 'ks', 'markerFaceColor', 'k')
plot([mean(includedAxonsThresh) + std(includedAxonsThresh), mean(includedAxonsThresh) - std(includedAxonsThresh)], [-2 -2], 'r', 'LineWidth', 2)
plot(mean(includedAxonsThresh), -2, 'rs', 'markerFaceColor', 'r')


hold off


%%
if 1    
    somaNormSD = [];
    axonNormSD = [];
    
    % response rate / curve plot
    figure('position', [100 100 650 400], 'color', [1 1 1])
    hold on
    for ii = 1:length(allSomasParams)
        if allSomasMaxResp(ii) >= respProbThresh
            %if max(cellInfo(ii).data(2,:)) > 0.5 % analyzed through at least prob = 0.5
            %xProj = cellInfo(ii).data(1,1):0.01:cellInfo(ii).data(1,end);
            xProj = 0:0.01:2;
            params = allSomasParams(ii,:);
            projection = 0.5 + 0.5*erf(params(1)*xProj + params(2));
            
            
            %plot(cellInfo(ii).data(1,:)/cellInfo(ii).thresh, cellInfo(ii).data(2,:), '.', 'color', midgetColor)
            plot(xProj/(-params(2)/params(1)), projection, 'Color', [0 0 0]);
            %         else
            %             disp(['excluded neuron ' num2str(cellInfo(ii).id)...
            %                 ' from response curve slope plots because it doesn''t reach 0.5 probability within analyzed region'])
            %         end
            
            somaNormSD = [somaNormSD -1/(sqrt(2)*params(2))];
            
        else
            disp(['excluded a soma because it didn''t reach ' num2str(respProbThresh) ' probability within analyzed region'])
        end
    end
    
    for ii = 1:length(allAxonsParams)
        if allAxonsMaxResp(ii) >= respProbThresh
            
            xProj = 0:0.01:2;
            params = allAxonsParams(ii,:);
            projection = 0.5 + 0.5*erf(params(1)*xProj + params(2));
            
            
            %plot(cellInfo(ii).data(1,:)/cellInfo(ii).thresh, cellInfo(ii).data(2,:), '.', 'color', midgetColor)
            plot(xProj/(-params(2)/params(1)), projection, 'Color', [1 0 0]);
            
            axonNormSD = [axonNormSD];
            axonNormSD = [axonNormSD -1/(sqrt(2)*params(2))];
                        
        else
            disp(['excluded an axon because it didn''t reach ' num2str(respProbThresh) ' probability within analyzed region'])
        end
    end
    
    hold off
    xlabel('normalized current amplitude', 'fontsize', 20)
    ylabel('response probability', 'fontsize', 20)
    set(gca, 'XLim', [0 2])
    set(gca, 'fontsize', 18)
    
    
%% 
    
    %histograms of normalized standard deviations
    binEdges = 0:0.05:0.5;

    figure('position', [300 300 300 300]); hA = axes; hold on
    psthPlotterBase(hA, binEdges', axonNormSD, 'lineColor', [1 0 0], 'fillHist', false)
    psthPlotterBase(hA, binEdges', somaNormSD, 'lineColor', [0 0 0], 'fillHist', false)
    xlabel('normalized SD')
    ylabel('number of cells')

end

%%

if 0
    
    absFig = figure;
    
    hold on
    
    plot(midgetSomaThresh, zeros(length(midgets), 1), 'ko')
    plot(midgetAxonThresh, zeros(length(midgets), 1), 'ro')
    
    plot(parasolSomaThresh, zeros(length(parasols), 1), 'k*')
    plot(parasolAxonThresh, zeros(length(parasols), 1), 'r*')
    
    plot(noAxonStimThresh, zeros(length(noAxonStimThresh), 1), 'ko')
    
    hold off
    
    relFig = figure;
    
    hold on
    
    plot(midgetAxonThresh./midgetSomaThresh, zeros(length(midgets), 1), 'ko')
    
    plot(parasolAxonThresh./parasolSomaThresh, zeros(length(parasols), 1), 'k*')
    
    hold off

end
