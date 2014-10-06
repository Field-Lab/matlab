
clear all
close all

%analysis parameters
recalcAll = false; %redo all response curve fits and standard deviation calculations
nBootstrapReps = 100;

analysisPathPrefix = '/snle/lab/Experiments/Array/Analysis/';

%arraySpacing = 60; %pitch of electrodes in microns

%colors!
blue = [50 70 247]/255; %blue
rust = [.8 .05 0.05];  %rust
grass = [90 156 0]/255; %pale grass
salmon = [255 124 59]/255;  %salmon

prefPW = 100;
%minPCrossing = 0.35; %minimum probability that must be crossed for cell to be included


%% cell/stim info

cellInfo = cell_list_sbc();

%%% separate out unstimulated cells into a list (remove from cellInfo struct)

count = 0;
for ii = length(cellInfo):-1:1
    if isempty(cellInfo(ii).stimElec)
        count = count + 1;
        cellInfoNoStim(count) = cellInfo(ii);
        cellInfo(ii) = [];
    end
end


%%


% %%% checking whether data is entered correctly
% 
% verMinCells = [];
% stimCells = [];
% for ii = 1:length(cellInfo)
%     if cellInfo(ii).verMin
%         verMinCells = [verMinCells cellInfo(ii).id];
%     end
%     if ~isempty(cellInfo(ii).stimElec)
%         stimCells = [stimCells cellInfo(ii).id];
%     end
% end


%% getting data from elecResp files

nCells = length(cellInfo);

for ii = 1:nCells
    if length(cellInfo(ii).PW) > 1
        cellInfo(ii).PW = prefPW;
%     else
%         PW = cellInfo(ii).PW;
    end
    if exist([analysisPathPrefix filesep cellInfo(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec) '_w' num2str(cellInfo(ii).PW) '.mat'], 'file')
        fullPath = [analysisPathPrefix filesep cellInfo(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec) '_w' num2str(cellInfo(ii).PW) '.mat'];
        load(fullPath)
    elseif exist([analysisPathPrefix filesep cellInfo(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec) '.mat'], 'file');
        fullPath = [analysisPathPrefix filesep cellInfo(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec) '.mat'];
        load(fullPath);
        disp(['no PW specifier in file elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec)])
    else
        error(['could not find file elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec)])
    end
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll);
    save(fullPath, 'elecResp')
    clear fullPath
        
    %get response rates
    cellInfo(ii).data = zeros(3, length(elecResp.stimInfo.stimAmps));
    cellInfo(ii).data(1,:) = abs(elecResp.stimInfo.stimAmps);
    cellInfo(ii).data(2,:) = elecResp.analysis.successRates;
    cellInfo(ii).data(3,:) = elecResp.stimInfo.nPulses;
    
    %remove unanalyzed data
    for jj = length(elecResp.stimInfo.stimAmps):-1:1
        if isempty(elecResp.analysis.type{jj})
            cellInfo(ii).data(:,jj) = [];
        end
    end
    
    %get fit parameters
    cellInfo(ii).params = elecResp.analysis.erfParams;
    cellInfo(ii).SD = 1/(sqrt(2)*cellInfo(ii).params(1)); %standard deviation of cumulative Gaussian fit, in µA
    cellInfo(ii).thresh = elecResp.analysis.threshold; % should equal -cellInfo(ii).params(2)/cellInfo(ii).params(1)
    
    %get ei information
    cellInfo(ii).ei = elecResp.cells.mainEI; %64 channels x 81 samples
    cellInfo(ii).eiAmps = max(-cellInfo(ii).ei, [], 2);
    
    %correct for change in gain between 2008 and 2010 datasets (from 840 to
    %440)  -- need to verify that this is the point in time when the switch
    %occurred!!
    if strcmp(cellInfo(ii).analysisPath(1:4), '2008')
        cellInfo(ii).ei = cellInfo(ii).ei*440/840;
        cellInfo(ii).eiAmps = cellInfo(ii).eiAmps*440/840;
    end
    
    clear elecResp
    
end

%% get information for unstimulated cells

for ii = 1:length(cellInfoNoStim)
    
    if length(cellInfoNoStim(ii).PW) > 1
        cellInfoNoStim(ii).PW = prefPW;
%     else
%         PW = cellInfoNoStim(ii).PW;
    end
    if exist([analysisPathPrefix filesep cellInfoNoStim(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfoNoStim(ii).id) '_p' num2str(cellInfoNoStim(ii).limStimElec) '_w' num2str(cellInfoNoStim(ii).PW) '.mat'], 'file')
        fullPath = [analysisPathPrefix filesep cellInfoNoStim(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfoNoStim(ii).id) '_p' num2str(cellInfoNoStim(ii).limStimElec) '_w' num2str(cellInfoNoStim(ii).PW) '.mat'];
        load(fullPath)
    elseif exist([analysisPathPrefix filesep cellInfoNoStim(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfoNoStim(ii).id) '_p' num2str(cellInfoNoStim(ii).limStimElec) '.mat'], 'file');
        fullPath = [analysisPathPrefix filesep cellInfoNoStim(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfoNoStim(ii).id) '_p' num2str(cellInfoNoStim(ii).limStimElec) '.mat'];
        load(fullPath);
        disp(['no PW specifier in file elecResp_n' num2str(cellInfoNoStim(ii).id) '_p' num2str(cellInfoNoStim(ii).limStimElec)])
    else
        error(['could not find file elecResp_n' num2str(cellInfoNoStim(ii).id) '_p' num2str(cellInfoNoStim(ii).limStimElec)])
    end
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll);
    save(fullPath, 'elecResp')
    clear fullPath
        
    
    %get response rates
    cellInfoNoStim(ii).data = zeros(3, length(elecResp.stimInfo.stimAmps));
    cellInfoNoStim(ii).data(1,:) = abs(elecResp.stimInfo.stimAmps);
    cellInfoNoStim(ii).data(2,:) = elecResp.analysis.successRates;
    cellInfoNoStim(ii).data(3,:) = elecResp.stimInfo.nPulses;
    
    %remove unanalyzed data
    for jj = length(elecResp.stimInfo.stimAmps):-1:1
        if isempty(elecResp.analysis.type{jj})
            cellInfoNoStim(ii).data(:,jj) = [];
        end
    end
    
    %determine minimum possible threshold
    minThreshData = cellInfoNoStim(ii).data;
    %add five data poinst with stim probability = 1 at 10% amplitude intervals
    for jj = 1:5
        x = size(minThreshData,2);
        minThreshData(1,x+1) = minThreshData(1,x)*1.1;
        minThreshData(2,x+1) = 1;
        minThreshData(3,x+1) = mean(cellInfoNoStim(ii).data(3,:));
        lockedAmps(x+1) = 0;
    end
    
    [params] = erfFitter(minThreshData, 2, -1, 'makePlot', 0, 'lockedAmps', lockedAmps);
    cellInfoNoStim(ii).minThresh = -params(2)/params(1);
    
    
    %get ei information
    cellInfoNoStim(ii).ei = elecResp.cells.mainEI; %64 channels x 81 samples
    cellInfoNoStim(ii).eiAmps = max(-cellInfoNoStim(ii).ei, [], 2);
    
    %eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile([analysisPathPrefix cellInfoNoStim(ii).eiPath]);
    %ei = eiFile.getImage(cellInfoNoStim(ii).id);
    %cellInfoNoStim(ii).ei = squeeze(ei(1,2:end,:));
    %cellInfoNoStim(ii).eiAmps = max(-cellInfoNoStim(ii).ei, [], 2);
    %eiFile.close()
    
    %correct for change in gain between 2008 and 2010 datasets (from 840 to
    %440) -- need to verify that this is the point in time when the switch
    %occurred!!
    if strcmp(cellInfoNoStim(ii).analysisPath(1:4), '2008')
        cellInfoNoStim(ii).ei = cellInfoNoStim(ii).ei*440/840;
        cellInfoNoStim(ii).eiAmps = cellInfoNoStim(ii).eiAmps*440/840;
    end
end


%% plot threshold vs. max negative peak ei signal

figure; hold on
for ii = 1:length(cellInfo)
    if cellInfo(ii).PW == 50
        if cellInfo(ii).verMin
            plot(max(cellInfo(ii).eiAmps), cellInfo(ii).thresh*cellInfo(ii).PW, 'o', 'MarkerEdgeColor', grass, 'MarkerFaceColor', grass, 'MarkerSize', 6)
        else
            plot(max(cellInfo(ii).eiAmps), cellInfo(ii).thresh*cellInfo(ii).PW, 'o', 'MarkerEdgeColor', grass, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
        end
    else %PW = 100
        if cellInfo(ii).verMin
            plot(max(cellInfo(ii).eiAmps), cellInfo(ii).thresh*cellInfo(ii).PW, 'o', 'MarkerEdgeColor', salmon, 'MarkerFaceColor', salmon, 'MarkerSize', 6)
        else
            plot(max(cellInfo(ii).eiAmps), cellInfo(ii).thresh*cellInfo(ii).PW, 'o', 'MarkerEdgeColor', salmon, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
        end
    end
end

for ii = 1:length(cellInfoNoStim)
    %plot(max(cellInfoNoStim(ii).eiAmps), 0, 'x', 'MarkerEdgeColor', 'black', 'MarkerSize', 6)

    if cellInfo(ii).PW == 50
        plot(max(cellInfoNoStim(ii).eiAmps)*[1 1], [cellInfoNoStim(ii).minThresh 10]*cellInfoNoStim(ii).PW,  'color', grass)
        plot(max(cellInfoNoStim(ii).eiAmps), cellInfoNoStim(ii).minThresh*cellInfoNoStim(ii).PW, 'x', 'MarkerEdgeColor', grass)
    else % PW = 100
        plot(max(cellInfoNoStim(ii).eiAmps)*[1 1], [cellInfoNoStim(ii).minThresh 10]*cellInfoNoStim(ii).PW,  'color', salmon)
        plot(max(cellInfoNoStim(ii).eiAmps), cellInfoNoStim(ii).minThresh*cellInfoNoStim(ii).PW, 'x', 'MarkerEdgeColor', salmon)
    end
end

text(110, 230, '50 µs PW', 'color', grass, 'fontweight', 'bold')
text(110, 215, '100 µs PW', 'color', salmon, 'fontweight', 'bold')


hold off
ylabel('threshold charge (pC)')
xlabel('peak negative ei signal (DAQ)')


set(gca, 'yLim', [0 50*ceil(0.02*max([cellInfo.thresh].*[cellInfo.PW]))])

