%% script to plot response curves of midgets and parasols from 2008-08-27-2-data002

clear all
close all

%analysis parameters
recalcAll = false; %redo all response curve fits and standard deviation calculations
nBootstrapReps = 100;
maxSuccRateCutoff = 0.4;


analysisPath = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/';
fileNameSuffix = ''; %use when _w50 or _w100 is included in filenames

pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';
datarunPath = '2008-08-27-2/data001-lh/data001-lh';
%exampleCells = [108 811 348 317]; %example curves in mosaic plot
examples.ids = [257 811 348 317]; %example curves in mosaic plot
examples.xLims = {[0 1.5],[0 1.5],[0 1.5],[0 1.5]};
%exampleXLims = {[0 1.5],[0 1.5],[0 1.5],[0 1.5]};

arraySpacing = 60; %pitch of electrodes in microns
cutoffMult = 0.5;

%colors!
blue = [50 70 247]/255; %blue
rust = [.8 .05 0.05];  %rust
grass = [90 156 0]/255; %pale grass
salmon = [255 124 59]/255;  %salmon

midgetColor = rust;
parasolColor = blue;
axonColor = grass;

cellTypes = {'onPar', 'onMidg'};

colors(1,:) = 1-(1-blue)*0.7; %ON parasol
colors(2,:) = 1-(1-rust)*0.7; %ON midget

%% cell/stim info

cellInfo = cell_list_2008_08_27_2();

%%% criteria for included cells: max (negative) signal must either
% (1) not fall on an edge electrode OR
% (2) be >= cutoffMult of the mean max signal of cells within class that meet (1)

[cellInfo cellInfoOffArray] = removeSmallSigEdgeCells(cellInfo, pathToEi, 'cutoffMult', 0.5);
%[cellInfo cellInfoOffArray] = removeSmallSigEdgeCells(cellInfo, pathToEi);

%%% checking whether data is entered correctly

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

%%% separate out unstimulated cells into a list (remove from cellInfo struct)

count = 0;
for ii = length(cellInfo):-1:1
    if isempty(cellInfo(ii).stimElec)
        count = count + 1;
        cellInfoNoStim(count) = cellInfo(ii);
        cellInfo(ii) = [];
    end
end


%% getting data from elecResp files

nCells = length(cellInfo);

for ii = 1:nCells
    load([analysisPath filesep 'elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec) fileNameSuffix '.mat'])
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll);
    save([analysisPath filesep 'elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec) fileNameSuffix '.mat'], 'elecResp')
    
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
end

%% apply maxSuccRateCutOff; add cells that don't pass to "no stim" category

count = 0;
for ii = length(cellInfo):-1:1
    if max(cellInfo(ii).data(2,:)) < maxSuccRateCutoff
        disp(['cell ' num2str(cellInfo(ii).id) ' didn''t meet success rate requirement of at least ' num2str(maxSuccRateCutoff)])
        
        %add cell to "cellInfoNoStim"
        ind = length(cellInfoNoStim)+1;
        cellInfoNoStim(ind).id = cellInfo(ii).id;
        cellInfoNoStim(ind).stimElec = cellInfo(ii).stimElec;
        cellInfoNoStim(ind).verMin = cellInfo(ii).verMin;
        cellInfoNoStim(ind).type = cellInfo(ii).type;
        
        %remove cell from "cellInfoNoStim"
        cellInfo(ii) = [];
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%% summary plot 

% load RFs
datarun = load_data(datarunPath);
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));
datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
datarun = get_sta_fits_from_vision(datarun,'all');

%%
params.threshBarLim = [0 250];
params.mosaicPlotLimX = [8.3 22.2];
params.mosaicPlotLimY = [9.6 23.45];
params.PW = 100;
params.chargeLim = 1.1043*params.PW;
params.safetyLim = 158.1;


params.yBarPos = [0.55 0.2];
params.OS = 0.12;

%cellTypeThreshCompareSummaryPlot(cellInfoPassesCutoff, cellInfoNoStim, datarun, colors, examples, params, cellTypes)

figure('position', [100 100 800 220], 'color', [1 1 1])
axes('position', [0.1 0.1 0.8 0.8])
make_1D_thresh_plot(cellInfo, colors, params, cellTypes)




%% response curve slope comparison plots

if 0
    % data for histogram plot
    midgAbove50normSD = []; %list standard deviation/threshold for each midget that is analyzable through prob = 0.5
    parAbove50normSD  = []; %list standard deviation/threshold for each midget that is analyzable through prob = 0.5
    
    % response rate / curve plot
    figure('position', [100 100 650 400], 'color', [1 1 1])
    hold on
    for ii = 1:nCells
        if max(cellInfo(ii).data(2,:)) > 0.5 % analyzed through at least prob = 0.5
            xProj = cellInfo(ii).data(1,1):0.01:cellInfo(ii).data(1,end);
            projection = 0.5 + 0.5*erf(cellInfo(ii).params(1)*xProj + cellInfo(ii).params(2));
            
            if strcmpi(cellInfo(ii).type, 'onMidg')
                plot(cellInfo(ii).data(1,:)/cellInfo(ii).thresh, cellInfo(ii).data(2,:), '.', 'color', midgetColor)
                plot(xProj/cellInfo(ii).thresh, projection, 'Color', midgetColor);
                midgAbove50normSD = [midgAbove50normSD cellInfo(ii).SD/cellInfo(ii).thresh];
            elseif strcmpi(cellInfo(ii).type, 'onPar')
                plot(cellInfo(ii).data(1,:)/cellInfo(ii).thresh, cellInfo(ii).data(2,:), '.', 'color', parasolColor)
                plot(xProj/cellInfo(ii).thresh, projection, 'Color', parasolColor);
                parAbove50normSD = [parAbove50normSD cellInfo(ii).SD/cellInfo(ii).thresh];
            else
                error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
            end
            
        else
            disp(['excluded neuron ' num2str(cellInfo(ii).id)...
                ' from response curve slope plots because it doesn''t reach 0.5 probability within analyzed region'])
        end
    end
    
    hold off
    xlabel('current amplitude, normalized to threshold', 'fontsize', 20)
    ylabel('response rate (%)', 'fontsize', 20)
    set(gca, 'XLim', [0 2])
    set(gca, 'fontsize', 18)
    
    
    % bar plot of standard deviations by cell type
    
    figure('position', [300 300 300 300])
    hold on
    plot(1, mean(midgAbove50normSD), 'o', 'MarkerFaceColor', midgetColor, 'MarkerEdgeColor', midgetColor)
    plot(2, mean(parAbove50normSD), 'o', 'MarkerFaceColor', parasolColor, 'MarkerEdgeColor', parasolColor)
    plot([1 1], [mean(midgAbove50normSD) - std(midgAbove50normSD), mean(midgAbove50normSD) + std(midgAbove50normSD)], 'Color', midgetColor)
    plot([2 2], [mean(parAbove50normSD) - std(parAbove50normSD), mean(parAbove50normSD) + std(parAbove50normSD)], 'Color', parasolColor)
    set(gca, 'xlim', [0 3], 'xTick', [1 2], 'xticklabel', {'ON midget', 'ON parasol'}, 'ylim', [0 0.5])
    ylabel('SD/mean of cumulative Gaussian fit')
    title(['normalized response curve slopes' 10 '(cells that can be analyzed at least to prob = 0.5)'])
    
    
    % histograms of normalized standard deviations
    histVector = 0:0.05:0.5;
    midgHist = hist(midgAbove50normSD, histVector);
    parHist = hist(parAbove50normSD, histVector);
    
    
    %grouped histogram
    figure('position', [300 300 300 300])
    histVals =  [midgHist; parHist];
    bar(histVector, histVals', 'grouped')
    set(gca, 'xlim', [0 0.5])
    xlabel('SD/mean of cumulative Gaussian fit')
    ylabel('number of cells')
    
    
    %overlapped (need to add transparencies)
    figure('position', [300 300 300 300]); hold on
    parH = bar(histVector, parHist');
    midgH = bar(histVector, midgHist');
    set(gca, 'xlim', [0 0.5])
    set(parH, 'faceColor', parasolColor, 'edgeColor', 'none')
    set(midgH, 'faceColor', midgetColor, 'edgeColor', 'none')
    xlabel('SD/mean of cumulative Gaussian fit')
    ylabel('number of cells')
    
end


%% plots of different variable-property relations

%% thresh vs. max spike amplitude

if 0
    
    eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
    for ii = 1:nCells
        ei = eiFile.getImage(cellInfo(ii).id);
        %calculate negative peak (absolute) signal on each electrode
        eiAmps = zeros(64, 1);
        for j = 1:64
            if ~(j==9||j==25||j==57)
                eiAmps(j) = max(max(-ei(1,j+1,:)));
            end
        end
        cellInfo(ii).maxSigAllElecs = max(eiAmps);
        cellInfo(ii).maxSigStimElec = eiAmps(cellInfo(ii).stimElec);
    end
    eiFile.close()
    
    %%% max amplitude on any electrode %%%
    figure; hold on
    for ii = 1:nCells
        if strcmpi(cellInfo(ii).type, 'onMidg')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', midgetColor, 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', parasolColor, 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        else
            error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
        end
    end
    xlabel('peak negative signal on any electrode')
    ylabel('threshold')
    hold off
    
    %%% max amplitude on any electrode %%%
    figure; hold on
    for ii = 1:nCells
        if strcmpi(cellInfo(ii).type, 'onMidg')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', midgetColor, 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', parasolColor, 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        else
            error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
        end
    end
    xlabel('peak negative signal on stimulation electrode')
    ylabel('threshold')
    hold off

end

%% thresh vs. max spike amplitude (stim elec)


%% thresh vs. COM and shifted COM

if 0
    
    checkForAxonSig = false;
    plotCOM = false;
    sigThresh = 0.25;
    
    % calculates COM
    eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
    [xCoords yCoords] = getElectrodeCoords61();
    
    for ii = 1:nCells
        ei = eiFile.getImage(cellInfo(ii).id);
        
        %calculate peak (absolute) signal on each electrode and normalize
        eiAmps = zeros(64, 1);
        for j = 1:64
            if ~(j==9||j==25||j==57)
                eiAmps(j) = max(max(abs(ei(1,j+1,:))));
            end
        end
        ei = ei/max(eiAmps);
        eiAmps = eiAmps/max(eiAmps);
        
        
        % cleaning out above-threshold axonal signals
        for jj = 1:length(cellInfo(ii).axonElecs)
            eiAmps(cellInfo(ii).axonElecs(jj)) = 0;
        end
        
        cellInfo(ii).COM = getCOM(ei, eiAmps, sigThresh, plotCOM, checkForAxonSig, cellInfo(ii).id); %in coordinates defined by getElectrodeCoords61()
        if checkForAxonSig
            uiwait()
        end
        
        % calculating distance from COM to stim electrode
        cellInfo(ii).stimDistToCOM = (arraySpacing/2)*norm(cellInfo(ii).COM - [xCoords(cellInfo(ii).stimElec); yCoords(cellInfo(ii).stimElec)]); %in microns, for 60-micron arrays
    end
    
    eiFile.close()
    
    
    
    %%% plot distance btwn stim electrode and COM vs. thresh %%%
    
    figure('position', [100 100 650 400], 'color', [1 1 1])
    
    hold on
    for ii = 1:nCells
        if strcmpi(cellInfo(ii).type, 'onMidg')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', midgetColor, 'MarkerSize', 6)
            else
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', parasolColor, 'MarkerSize', 6)
            else
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        else
            error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
        end
    end
    hold off
    xlabel(['distance between stimulating electrode and', 10, 'electrical image center of mass (\mum)'], 'fontsize', 20)
    ylabel('threshold (uA)', 'fontsize', 20)
    set(gca, 'xlim', [-10 100], 'ylim', [0 2], 'xtick', [0 20 40 60 80 100], 'fontsize', 18)
    
    
    
    %%% plot thresh vs distance btwn stim elec and different offsets from COM in axon direction %%%
    
    %approximated angle
    axonAngle = 300; %in getElectrodeCoords61() coordinates: elec 13 on +x axis and elec 62 on +y axis
    
    shiftLengths = [0 10 20 30 40 50 60]; %distance from COM, in microns
    
    distsToOffsetCOM = zeros(nCells, length(shiftLengths));
    for ii = 1:nCells
        for jj = 1:length(shiftLengths)
            offsetCOM = cellInfo(ii).COM*arraySpacing/2 + shiftLengths(jj)*[cos(axonAngle*pi/180); sin(axonAngle*pi/180)]; %in microns
            distsToOffsetCOM(ii,jj) = norm(offsetCOM - (arraySpacing/2)*[xCoords(cellInfo(ii).stimElec); yCoords(cellInfo(ii).stimElec)]);
        end
    end
    
    for jj = 1:length(shiftLengths)
        figure
        hold on
        for ii = 1:nCells
            if strcmpi(cellInfo(ii).type, 'onMidg')
                if cellInfo(ii).verMin
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', midgetColor, 'MarkerSize', 6)
                else
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', midgetColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
                end
            elseif strcmpi(cellInfo(ii).type, 'onPar')
                if cellInfo(ii).verMin
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', parasolColor, 'MarkerSize', 6)
                else
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', parasolColor, 'MarkerFaceColor', 'none', 'MarkerSize', 6)
                end
            else
                error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
            end
        end
        
        xlabel('distance between stim elec and shifted COM')
        ylabel('threshold')
        title(['shift length = ' num2str(shiftLengths(jj)) 'µm (in direction of axon)'])
    end
    
end

%% old figures

% % midgets
% midgetColors = hsv(nMidgets);
% 
% figure
% subplot(2,1,1)
% hold on
% for i = 1:nMidgets
%    plot(midgetStimAmps{i}, midgetResponses{i}*100, '.', 'MarkerEdgeColor', midgetColors(i,:),...
%        'MarkerFaceColor', midgetColors(i,:), 'MarkerSize', 10);
%    current = plot(midgetProjections{i}(1,:), midgetProjections{i}(2,:)*100);
%    set(findobj(current,'Type','line'),'Color', midgetColors(i,:))
% end
% hold off
% title('On Midget Response Curves')
% xlabel('stimulation current amplitude (\muA)')
% ylabel('percent of pulses that elicit a response')
% 
% 
% 
% % parasols
% parasolColors = hsv(nParasols);
% 
% subplot(2,1,2)
% hold on
% for i = 1:nParasols
%     plot(parasolStimAmps{i}, parasolResponses{i}*100, '.', 'MarkerEdgeColor', parasolColors(i,:),...
%         'MarkerFaceColor', parasolColors(i,:), 'MarkerSize', 10);
%     current = plot(parasolProjections{i}(1,:), parasolProjections{i}(2,:)*100);
%     set(findobj(current,'Type','line'),'Color', parasolColors(i,:))
% end
% hold off
% title('On Parasol Response Curves')
% xlabel('stimulation current amplitude (\muA)')
% ylabel('percent of pulses that elicit a response')
