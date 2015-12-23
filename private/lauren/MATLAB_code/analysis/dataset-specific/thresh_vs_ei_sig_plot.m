%parameters for all datasets

clear all

%analysis parameters
recalcAll = false; %redo all response curve fits and standard deviation calculations
nBootstrapReps = 100;
cutoffMult = 0.5;
maxSuccRateCutoff = 0.4;
exclude30ArrayData = true;

cellTypes = {'onPar', 'offPar', 'onMidg', 'offMidg', 'sbc'};

%colors!
blue = [50 70 247]/255; %blue
rust = [.8 .05 0.05];  %rust
grass = [90 156 0]/255; %pale grass
salmon = [255 124 59]/255;  %salmon

colors(1,:) = 1-(1-blue)*0.7; %ON parasol
colors(2,:) = blue*0.7; %OFF parasol
colors(3,:) = 1-(1-rust)*0.7; %ON midget
colors(4,:) = rust*0.7; %OFF midget
colors(5,:) = grass; %SBC

elecAreaDS4.means = [];
elecAreaDS4.analysisPaths = {};

%% dataset details

analysisPath{1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/';
fileNameSuffix{1} = ''; %use when _w50 or _w100 is included in filenames
PW{1} = 100;
pathToEi{1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';
cellInfo{1} = cell_list_2008_08_27_2();

analysisPath{2} = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/';
fileNameSuffix{2} = '_w50'; %use when _w50 or _w100 is included in filenames
PW{2} = 50;
pathToEi{2} = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030-lh/data030-lh.ei';
cellInfo{2} = cell_list_2011_01_11_0();

analysisPath{3} = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data002/';
fileNameSuffix{3} = '_w50'; %use when _w50 or _w100 is included in filenames
PW{3} = 50;
pathToEi{3} = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data000/data000.ei';
cellInfo{3} = cell_list_2011_10_25_4();
%add fields to make consistent with DS4
for ii = 1:length(cellInfo{3})
    cellInfo{3}(ii).analysisPath = '2011-10-25-4/data002';
    cellInfo{3}(ii).meanElecArea = 128.7;
end

%this one's a little different...collection of SBCs from various preps
analysisPath{4} = '/snle/lab/Experiments/Array/Analysis/';
prefPW = 100;
cellInfo{4} = cell_list_sbc(exclude30ArrayData);

% remove cells that are in dataset 2011-10-25-4/data002 because they're
% already accounted for in dataset #3
for ii = length(cellInfo{4}):-1:1
    if strcmpi(cellInfo{4}(ii).analysisPath, '2011-10-25-4/data002')
        disp(['removed cell ' num2str(cellInfo{4}(ii).id) ' from fourth dataset because it''s already in 3rd dataset'])
        cellInfo{4}(ii) = [];
    end
end

%% apply exclusions

%%% separate out off-array cells for datasets corresponding to single preparation
for iDS = 1:3
    [cellInfo{iDS} cellInfoOffArray{iDS}] = removeSmallSigEdgeCells(cellInfo{iDS}, pathToEi{iDS}, 'cutOffMult', cutoffMult);
end

for iDS = 1:4
    %make sure all remaining cells have been analyzed
    for ii = 1:length(cellInfo{iDS})
        if cellInfo{iDS}(ii).stimElec == 0 %indicates unfinished analysis
            error(['unfinished analysis for cell ' num2str(cellInfo{iDS}(ii).id) ' from ' cellInfo{iDS}(ii).analysisPath])
        end
    end
    
    %%% separate out unstimulated cells
    count = 0;
    for ii = length(cellInfo{iDS}):-1:1
        if isempty(cellInfo{iDS}(ii).stimElec)
            count = count + 1;
            cellInfoNoStim{iDS}(count) = cellInfo{iDS}(ii);
            cellInfo{iDS}(ii) = [];
        end
    end
end

%% get analysis

for iDS = 1:3
    % load analysis from elecResp files
    for ii = 1:length(cellInfo{iDS})
        load([analysisPath{iDS} filesep 'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) fileNameSuffix{iDS} '.mat'])
        elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll);
        save([analysisPath{iDS} filesep 'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) fileNameSuffix{iDS} '.mat'], 'elecResp')
        
        cellInfo{iDS}(ii).params = elecResp.analysis.erfParams;
        cellInfo{iDS}(ii).SD = 1/(sqrt(2)*cellInfo{iDS}(ii).params(1)); %standard deviation of cumulative Gaussian fit, in µA
        cellInfo{iDS}(ii).thresh = elecResp.analysis.threshold; % should equal -cellInfo{iDS}(ii).params(2)/cellInfo{iDS}(ii).params(1)
        cellInfo{iDS}(ii).chargeThresh = PW{iDS}*cellInfo{iDS}(ii).thresh;
        
        
        %get maximum measured response rate
        succRates = elecResp.analysis.successRates;
        %remove unanalyzed data
        for jj = length(elecResp.stimInfo.stimAmps):-1:1
            if isempty(elecResp.analysis.type{jj})
                succRates(jj) = [];
            end
        end
        cellInfo{iDS}(ii).maxSuccRate = max(succRates);
        
        clear elecResp succRates
    end
    
    %correct for change in gain
    for ii = 1:length(cellInfo{iDS})
        cellInfo{iDS}(ii).maxSigAllElecs = 440/840*cellInfo{iDS}(ii).maxSigAllElecs;
    end
    
    for ii = 1:length(cellInfoNoStim{iDS})
        cellInfoNoStim{iDS}(ii).maxSigAllElecs = 440/840*cellInfoNoStim{iDS}(ii).maxSigAllElecs;
    end
end


iDS = 4;
for ii = 1:length(cellInfo{iDS})
    if length(cellInfo{iDS}(ii).PW) > 1
        
        %use other PW if preferred PW doesn't fulfill maxSuccRateCutoff
        tmp = load([analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '_w' num2str(prefPW) '.mat']);
        
        %determine maximum measured response rate
        succRatesTmp = tmp.elecResp.analysis.successRates;
        for jj = length(tmp.elecResp.stimInfo.stimAmps):-1:1
            if isempty(tmp.elecResp.analysis.type{jj})
                succRatesTmp(jj) = [];
            end
        end
        
        if max(succRatesTmp) >= maxSuccRateCutoff
            cellInfo{iDS}(ii).PW = prefPW;
        else
            cellInfo{iDS}(ii).PW = cellInfo{iDS}(ii).PW(cellInfo{iDS}(ii).PW~=prefPW);
            disp(['Preferred pulse width doesn''t satisfy maxSuccRateCutoff criterion for cell ' num2str(cellInfo{iDS}(ii).id) ': using other pulse width (' num2str(cellInfo{iDS}(ii).PW) ')'])
        end
    end
    if exist([analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '_w' num2str(cellInfo{iDS}(ii).PW) '.mat'], 'file')
        fullPath = [analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '_w' num2str(cellInfo{iDS}(ii).PW) '.mat'];
        load(fullPath)
    elseif exist([analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '.mat'], 'file');
        fullPath = [analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
            'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '.mat'];
        load(fullPath);
        disp(['no PW specifier in file elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec)])
    else
        error(['could not find file elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec)])
    end
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll);
    save(fullPath, 'elecResp')
    clear fullPath
    
    %get fit parameters
    cellInfo{iDS}(ii).params = elecResp.analysis.erfParams;
    cellInfo{iDS}(ii).SD = 1/(sqrt(2)*cellInfo{iDS}(ii).params(1)); %standard deviation of cumulative Gaussian fit, in µA
    cellInfo{iDS}(ii).thresh = elecResp.analysis.threshold; % should equal -cellInfo{iDS}(ii).params(2)/cellInfo{iDS}(ii).params(1)
    cellInfo{iDS}(ii).chargeThresh = cellInfo{iDS}(ii).PW*cellInfo{iDS}(ii).thresh;

    %get maximum measured response rate
    succRates = elecResp.analysis.successRates;
    %remove unanalyzed data
    for jj = length(elecResp.stimInfo.stimAmps):-1:1
        if isempty(elecResp.analysis.type{jj})
            succRates(jj) = [];
        end
    end
    cellInfo{iDS}(ii).maxSuccRate = max(succRates);
    
    %get ei information
    cellInfo{iDS}(ii).ei = elecResp.cells.mainEI; %64 channels x 81 samples
    cellInfo{iDS}(ii).eiAmps = max(-cellInfo{iDS}(ii).ei, [], 2);
    
    %correct for change in gain between 2008 and 2010 datasets (from 840 to
    %440)  -- need to verify that all 2008 datasets are gain = 440
    if strcmp(cellInfo{iDS}(ii).analysisPath(1:4), '2008')
        cellInfo{iDS}(ii).ei = cellInfo{iDS}(ii).ei*440/840;
        cellInfo{iDS}(ii).eiAmps = cellInfo{iDS}(ii).eiAmps*440/840;
    end
    
    cellInfo{iDS}(ii).maxSigAllElecs = max(cellInfo{iDS}(ii).eiAmps);
    
    clear elecResp succRate
end


%% apply maxSuccRateCutOff

for jj = 1:length(cellInfo)
    count = 0;
    for ii = length(cellInfo{jj}):-1:1
        if cellInfo{jj}(ii).maxSuccRate < maxSuccRateCutoff
            count = count + 1;
            cellInfoFailsCutoff{jj}(count) = cellInfo{jj}(ii);
            cellInfo{jj}(ii) = [];
            disp(['cell ' num2str(cellInfoFailsCutoff{jj}(count).id) ' didn''t meet success rate requirement of at least ' num2str(maxSuccRateCutoff)])
        end
    end
end


%% plot all sbc thresholds in same style as within-piece comparisons
figure('position', [100 100 800 220], 'color', [1 1 1])
axes('position', [0.1 0.1 0.8 0.8])

params.threshBarLim = [0 250];
params.PW = [50 100];
params.yBarPos = [0.2];
params.OS = 0.12;

includeDS3 = true; %whether or not to include data from 2011-10-25-4, which are also also plotted in the within-piece comparison
include_2012_01_27_3 = false; %alternative within-piece comparison data

%collect SBC info for 1D plot
count = 0;
if includeDS3
    for ii=1:length(cellInfo{3})
        if strcmpi(cellInfo{3}(ii).type, 'sbc')
            count = count+1;
            cellInfoFor1DPlot(count).id = cellInfo{3}(ii).id;
            cellInfoFor1DPlot(count).type = 'sbc';
            cellInfoFor1DPlot(count).verMin = cellInfo{3}(ii).verMin;
            cellInfoFor1DPlot(count).thresh = cellInfo{3}(ii).thresh;
            cellInfoFor1DPlot(count).PW = PW{3};
            cellInfoFor1DPlot(count).analysisPath = cellInfo{3}(ii).analysisPath;
            cellInfoFor1DPlot(count).meanElecArea = cellInfo{3}(ii).meanElecArea;
        end
    end
end
for ii = 1:length(cellInfo{4})
    
    if include_2012_01_27_3 || ~strcmpi(cellInfo{4}(ii).analysisPath, '2012-01-27-3/data002')
        count = count+1;
        cellInfoFor1DPlot(count).id = cellInfo{4}(ii).id;
        cellInfoFor1DPlot(count).type = 'sbc';
        cellInfoFor1DPlot(count).verMin = cellInfo{4}(ii).verMin;
        cellInfoFor1DPlot(count).thresh = cellInfo{4}(ii).thresh;
        cellInfoFor1DPlot(count).PW = cellInfo{4}(ii).PW;
        cellInfoFor1DPlot(count).analysisPath = cellInfo{4}(ii).analysisPath;
        cellInfoFor1DPlot(count).meanElecArea = cellInfo{4}(ii).meanElecArea;
    end
end

%calculate mean mean electrode area across arrays used
pieces = {};
meanElecAreas = [];
for ii = 1:length(cellInfoFor1DPlot)
    pieceNo = cellInfoFor1DPlot(ii).analysisPath(1:strfind(cellInfoFor1DPlot(ii).analysisPath,filesep)-1);
    if sum(cellfun(@(s)(strcmp(s, pieceNo)), pieces))==0 %haven't come across this piece yet
        pieces{end+1} = pieceNo;
        
        %make sure mean electrode area has been specified
        if ~isempty(cellInfoFor1DPlot(ii).meanElecArea)
            meanElecAreas = [meanElecAreas cellInfoFor1DPlot(ii).meanElecArea];
        else
            error('the mean electrode area has not been specified for all cells being plotted')
        end
    end
end

params.safetyLim = mean(meanElecAreas); %in pC

make_1D_thresh_plot(cellInfoFor1DPlot, colors(5,:), params, {'sbc'})


%% gather number of unstimulated or maxSuccRateCutOff failures for each
% cell type (for calculation of quartiles)

nNoStimByType = zeros(1,length(cellTypes));

for kk = 1:length(cellInfoFailsCutoff)
    for ii = 1:length(cellInfoFailsCutoff{kk})
        foundMatch = false;
        for jj = 1:length(cellTypes)
            if strcmpi(cellTypes{jj}, cellInfoFailsCutoff{kk}(ii).type)
                nNoStimByType(jj) = nNoStimByType(jj) + 1;
                foundMatch = true;
                break
            end
        end
        if ~foundMatch
            error(['unrecognized cell type for neuron ' num2str(cellInfoFailsCutoff{kk}(ii).id)])
        end
    end
end

for kk = 1:length(cellInfoNoStim)
    for ii = 1:length(cellInfoNoStim{kk})
        foundMatch = false;
        for jj = 1:length(cellTypes)
            if strcmpi(cellTypes{jj}, cellInfoNoStim{kk}(ii).type)
                nNoStimByType(jj) = nNoStimByType(jj) + 1;
                foundMatch = true;
                break
            end
        end
        if ~foundMatch
            error(['unrecognized cell type for neuron ' num2str(cellInfoNoStim{kk}(ii).id)])
        end
    end
end


%% plot thresh vs. max amplitude on any electrode %% (and gather thresholds by cell type)

for ii = 1:length(cellTypes)
    thresh{ii}.ver = zeros(0,2);
    thresh{ii}.unver = zeros(0,2);
end

figure; hold on

for kk = 1:length(cellInfo)
    for ii = 1:length(cellInfo{kk})
        foundMatch = false;
        for jj = 1:length(cellTypes)
            
            if strcmpi(cellTypes{jj}, cellInfo{kk}(ii).type)
                if cellInfo{kk}(ii).verMin
                    plot(cellInfo{kk}(ii).maxSigAllElecs, cellInfo{kk}(ii).chargeThresh, 'o', 'MarkerEdgeColor', colors(jj,:), 'MarkerFaceColor', colors(jj,:), 'MarkerSize', 6)
                    thresh{jj}.ver = [thresh{jj}.ver; zeros(1,2)];
                    thresh{jj}.ver(end,1) = cellInfo{kk}(ii).maxSigAllElecs;
                    thresh{jj}.ver(end,2) = cellInfo{kk}(ii).chargeThresh;
                else
                    plot(cellInfo{kk}(ii).maxSigAllElecs, cellInfo{kk}(ii).chargeThresh, 'o', 'MarkerEdgeColor', colors(jj,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
                    thresh{jj}.unver = [thresh{jj}.unver; zeros(1,2)];
                    thresh{jj}.unver(end,1) = cellInfo{kk}(ii).maxSigAllElecs;
                    thresh{jj}.unver(end,2) = cellInfo{kk}(ii).chargeThresh;
                end
                
                foundMatch = true;
                break
            end
        end
        
        if ~foundMatch
            error(['unrecognized cell type for neuron ' num2str(cellInfo{kk}(ii).id)])
        end
    end
end

xlabel('peak negative signal on any electrode')
ylabel('threshold')

%% plot thresh summaries by type
for ii = 1:length(cellTypes)
    threshAll = [thresh{ii}.unver(:,2); thresh{ii}.ver(:,2)];
    threshAllForQuartiles = [threshAll; 10*max(threshAll)*ones(nNoStimByType(ii),1)]; %set the unknown values to 10x the max known value
    
    quartiles = quantile(threshAllForQuartiles, [0.25 0.5 0.75]);
    
    if quartiles(2) > max(threshAll)
        error(['Not enough known data to compute median for cell type ' cellTypes{ii}])
    else
        plot(400+30*ii, quartiles(2) , 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:))
    end
    if quartiles(3) > max(threshAll)
        warning(['Not enough known data to compute third quartile for cell type ' cellTypes{ii}])
        plot((400+30*ii)*[1 1], [quartiles(1) quartiles(2)], '-', 'Color', colors(ii,:))
    else
        plot((400+30*ii)*[1 1], [quartiles(1) quartiles(3)], '-', 'Color', colors(ii,:))
    end
    
    %plot mean +/- SD of 'verified' group (excludes unstimulated cells)
    %plot((400+30*ii)*[1 1], mean(thresh{ii}.ver(:,2)) + [-std(thresh{ii}.ver(:,2)) std(thresh{ii}.ver(:,2))], '-', 'Color', colors(ii,:))
    %plot(400+30*ii, mean(thresh{ii}.ver(:,2)), 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:))
    
    %threshAll = [thresh{ii}.unver(:,2); thresh{ii}.ver(:,2)];
    %plot((410+30*ii)*[1 1], mean(threshAll) + [-std(threshAll) std(threshAll)], '-', 'Color', colors(ii,:))
    %plot(410+30*ii, mean(threshAll), 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', [1 1 1] - 0.5*([1 1 1]-colors(ii,:)))
end


hold off


%% plot thresh summaries by type in same style as within-piece comparisons
figure('position', [100 100 800 220], 'color', [1 1 1])
axes('position', [0.1 0.1 0.8 0.8]); hold on

threshBarLim = [0 250];
yBarPos = [1.25 1 0.75 0.5 0.25];

plot([threshBarLim(1) threshBarLim(2)], [0 0], 'k-', 'Linewidth', 1)
plot([threshBarLim(1) threshBarLim(1)], [-0.1 0.1], 'k-', 'Linewidth', 1)
plot([threshBarLim(2) threshBarLim(2)], [-0.1 0.1], 'k-', 'Linewidth', 1)


for ii = 1:length(cellTypes)
    threshAll = [thresh{ii}.unver(:,2); thresh{ii}.ver(:,2)];
    threshAllForQuartiles = [threshAll; 10*max(threshAll)*ones(nNoStimByType(ii),1)]; %set the unknown values to 10x the max known value
    
    quartiles = quantile(threshAllForQuartiles, [0.25 0.5 0.75]);
    
    if quartiles(2) > max(threshAll)
        error(['Not enough known data to compute median for cell type ' cellTypes{ii}])
    else
        plot(quartiles(2), yBarPos(ii) , 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:))
    end
    if quartiles(3) > max(threshAll)
        warning(['Not enough known data to compute third quartile for cell type ' cellTypes{ii}])
        plot([quartiles(1) quartiles(2)], yBarPos(ii)*[1 1], '-', 'Color', colors(ii,:))
    else
        plot([quartiles(1) quartiles(3)], yBarPos(ii)*[1 1], '-', 'Color', colors(ii,:))
    end
    
    %plot mean +/- SD of 'verified' group (excludes unstimulated cells)
    %plot((400+30*ii)*[1 1], mean(thresh{ii}.ver(:,2)) + [-std(thresh{ii}.ver(:,2)) std(thresh{ii}.ver(:,2))], '-', 'Color', colors(ii,:))
    %plot(400+30*ii, mean(thresh{ii}.ver(:,2)), 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:))
    
    %threshAll = [thresh{ii}.unver(:,2); thresh{ii}.ver(:,2)];
    %plot((410+30*ii)*[1 1], mean(threshAll) + [-std(threshAll) std(threshAll)], '-', 'Color', colors(ii,:))
    %plot(410+30*ii, mean(threshAll), 's', 'MarkerEdgeColor', colors(ii,:),
    %'MarkerFaceColor', [1 1 1] - 0.5*([1 1 1]-colors(ii,:)))
end

text(0, -0.5, '0 (0)', 'HorizontalAlignment', 'center', 'FontSize', 18)
text(0.5*threshBarLim(2), -0.5, [num2str(0.5*threshBarLim(2))], 'HorizontalAlignment', 'center', 'FontSize', 18);
text(threshBarLim(2), -0.5, [num2str(threshBarLim(2))], 'HorizontalAlignment', 'center', 'FontSize', 18);
text(0.5*threshBarLim(2), -1, 'threshold: pC', 'HorizontalAlignment', 'center', 'FontSize', 20);
hold off
set(gca, 'YLim', [-1 1.5], 'XLim', [threshBarLim(1) threshBarLim(2)])
axis off



%% max amplitude on stim electrode %%%

% for kk = 1:length(pathToEi)
% eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi{iDS});
% for ii = 1:nCells
%     ei = eiFile.getImage(cellInfo(ii).id);
%     %calculate negative peak (absolute) signal on each electrode
%     eiAmps = zeros(64, 1);
%     for j = 1:64
%         if ~(j==9||j==25||j==57)
%             eiAmps(j) = max(max(-ei(1,j+1,:)));
%         end
%     end
%     cellInfo(ii).maxSigAllElecs = max(eiAmps);
%     cellInfo(ii).maxSigStimElec = eiAmps(cellInfo(ii).stimElec);
% end
% eiFile.close()
% end


% figure; hold on
% for ii = 1:nCells
%     if strcmpi(cellInfo(ii).type, 'onMidg')
%         if cellInfo(ii).verMin
%             plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', midgetColor, 'MarkerSize', 6)
%         else
%             plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
%         end
%     elseif strcmpi(cellInfo(ii).type, 'onPar')
%         if cellInfo(ii).verMin
%             plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', parasolColor, 'MarkerSize', 6)
%         else
%             plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
%         end
%     else
%         error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
%     end
% end
% xlabel('peak negative signal on stimulation electrode')
% ylabel('threshold')
% hold off
