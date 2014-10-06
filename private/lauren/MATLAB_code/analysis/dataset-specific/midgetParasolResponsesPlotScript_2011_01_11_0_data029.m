%% script to plot response curves of midgets and parasols from 2011-01-11-0-data029

clear all
close all

%analysis parameters
recalcAll = false; %redo all response curve fits and standard deviation calculations
nBootstrapReps = 100;

maxSuccRateCutoff = 0.4;

analysisPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/';
fileNameSuffix = '_w50'; %use when _w50 or _w100 is included in filenames

pathToEi = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030-lh/data030-lh.ei';
datarunPath = '2011-01-11-0/data030-lh/data030-lh';

examples.ids = [407 617 302 680];
%examples.ids = [407 617 302 754];
%examples.ids = [152 888 767 338]; %example curves in mosaic plot: ON par, OFF par, ON midg, OFF midg
%alternatives: [407/2 392 631]

examples.xLims = {[0 2],[0 2],[0 2],[0 2]};

arraySpacing = 60; %pitch of electrodes in microns
cutoffMult = 0.5; %determines which cells will be included in analysis

cellTypes = {'onPar', 'offPar', 'onMidg', 'offMidg'};

%colors!
blue = [50 70 247]/255;
rust = [.8 .05 0.05];
grass = [90 156 0]/255;
salmon = [255 124 59]/255;

%midgetColor = rust;
%parasolColor = blue;
%axonColor = grass;

colors = zeros(4,3);

% colors(1,:) = 1-(1-blue)*0.7; %ON parasol
% colors(2,:) = blue*0.7; %OFF parasol
% colors(3,:) = 1-(1-rust)*0.7; %ON midget
% colors(4,:) = rust*0.7; %OFF midget


colors(1,:) = [0 0 0];
colors(2,:) = [0 1 1];
colors(3,:) = [0.2 1 0];
colors(4,:) = [1 0 0];



%% cell/stim info

cellInfo = cell_list_2011_01_11_0();


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


%%% criteria for included cells: max (negative) signal must either
% (1) not fall on an edge electrode OR
% (2) be >= cutoffMult of the mean max signal of cells within class that meet (1)

[cellInfo cellInfoOffArray] = removeSmallSigEdgeCells(cellInfo, pathToEi, 'cutoffMult', 0.5);
%[cellInfo cellInfoOffArray] = removeSmallSigEdgeCells(cellInfo, pathToEi);

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

%% apply maxSuccRateCutOff

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
datarun = load_index(datarun);
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));
datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
%datarun = get_sta_fits_from_vision(datarun,'all');

datarun.names.obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030-lh';
datarun = load_obvius_sta_fits(datarun);
datarun = get_sta_fits_from_obvius(datarun, 'all');

datarun = compute_monitor_to_array_transformation(datarun);



%%

params.threshBarLim = [0 250];
params.mosaicPlotLimX = [6 26];
params.mosaicPlotLimY = [6 26];
params.PW = 50;
params.yBarPos = [1.25 0.9 0.55 0.2];
params.OS = 0.12;
params.mosaicLineWidths = [2 2 2 2];
params.chargeLim = 4.3913*params.PW;
params.safetyLim = 166.3;


colors(1,:) = 1-(1-blue)*0.7; %ON parasol
colors(2,:) = blue*0.7; %OFF parasol
colors(3,:) = 1-(1-rust)*0.7; %ON midget
colors(4,:) = rust*0.7; %OFF midget

cellTypeThreshCompareSummaryPlot(cellInfo, cellInfoNoStim, datarun, colors, examples, params, cellTypes)

figure('position', [100 100 800 220], 'color', [1 1 1])
axes('position', [0.1 0.1 0.8 0.8])
make_1D_thresh_plot(cellInfo, colors, params, cellTypes, examples.ids)



%% response curve slope comparison plots

if 0
    
    % data for histogram plot
    onMidgAbove50normSD = []; %list standard deviation/threshold for each midget that is analyzable through prob = 0.5
    offMidgAbove50normSD = [];
    onParAbove50normSD  = [];
    offParAbove50normSD = [];

    
    % response rate / curve plot
    figure('position', [100 100 650 400], 'color', [1 1 1])
    hold on
    for ii = 1:nCells
        if max(cellInfo(ii).data(2,:)) > 0.5 % analyzed through at least prob = 0.5
            xProj = cellInfo(ii).data(1,1):0.01:cellInfo(ii).data(1,end);
            projection = 0.5 + 0.5*erf(cellInfo(ii).params(1)*xProj + cellInfo(ii).params(2));
            
            if strcmpi(cellInfo(ii).type, 'onMidg')
                plot(cellInfo(ii).data(1,:)/cellInfo(ii).thresh, cellInfo(ii).data(2,:), '.', 'color', colors(3,:))
                plot(xProj/cellInfo(ii).thresh, projection, 'Color', colors(3,:));
                onMidgAbove50normSD = [onMidgAbove50normSD cellInfo(ii).SD/cellInfo(ii).thresh]; %#ok<*AGROW>
            elseif strcmpi(cellInfo(ii).type, 'offMidg')
                plot(cellInfo(ii).data(1,:)/cellInfo(ii).thresh, cellInfo(ii).data(2,:), '.', 'color', colors(4,:))
                plot(xProj/cellInfo(ii).thresh, projection, 'Color', colors(4,:));
                offMidgAbove50normSD = [offMidgAbove50normSD cellInfo(ii).SD/cellInfo(ii).thresh];
            elseif strcmpi(cellInfo(ii).type, 'onPar')
                plot(cellInfo(ii).data(1,:)/cellInfo(ii).thresh, cellInfo(ii).data(2,:), '.', 'color', colors(1,:))
                plot(xProj/cellInfo(ii).thresh, projection, 'Color', colors(1,:));
                onParAbove50normSD = [onParAbove50normSD cellInfo(ii).SD/cellInfo(ii).thresh];
            elseif strcmpi(cellInfo(ii).type, 'offPar')
                plot(cellInfo(ii).data(1,:)/cellInfo(ii).thresh, cellInfo(ii).data(2,:), '.', 'color', colors(2,:))
                plot(xProj/cellInfo(ii).thresh, projection, 'Color', colors(2,:));
                offParAbove50normSD = [offParAbove50normSD cellInfo(ii).SD/cellInfo(ii).thresh];
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
    plot(1, mean(onMidgAbove50normSD), 'o', 'MarkerFaceColor', colors(3,:), 'MarkerEdgeColor', colors(3,:))
    plot(1.5, mean(offMidgAbove50normSD), 'o', 'MarkerFaceColor', colors(4,:), 'MarkerEdgeColor', colors(4,:))

    plot(2, mean(onParAbove50normSD), 'o', 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:))
    plot(2.5, mean(offParAbove50normSD), 'o', 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', colors(2,:))

    
    plot([1 1],     [mean(onMidgAbove50normSD)  - std(onMidgAbove50normSD),  mean(onMidgAbove50normSD)  + std(onMidgAbove50normSD)],  'Color', colors(3,:))
    plot([1.5 1.5], [mean(offMidgAbove50normSD) - std(offMidgAbove50normSD), mean(offMidgAbove50normSD) + std(offMidgAbove50normSD)], 'Color', colors(4,:))

    plot([2 2],     [mean(onParAbove50normSD)  - std(onParAbove50normSD),  mean(onParAbove50normSD)  + std(onParAbove50normSD)],  'Color', colors(1,:))
    plot([2.5 2.5], [mean(offParAbove50normSD) - std(offParAbove50normSD), mean(offParAbove50normSD) + std(offParAbove50normSD)],  'Color', colors(2,:))

    set(gca, 'xlim', [0.5 3], 'xTick', [1 1.5 2 2.5], 'xticklabel', {'ON-M', 'OFF-M', 'ON-P', 'OFF-P'}, 'ylim', [0 0.5])
    ylabel('SD/mean of cumulative Gaussian fit')
    title(['normalized response curve slopes' 10 '(cells that can be analyzed at least to prob = 0.5)'])
    
    
    % histograms of normalized standard deviations
    histVector = 0:0.05:0.5;
    onMidgHist =  hist(onMidgAbove50normSD, histVector);
    offMidgHist = hist(offMidgAbove50normSD, histVector);
    onParHist =  hist(onParAbove50normSD, histVector);
    offParHist = hist(offParAbove50normSD, histVector);
    
    %grouped histogram
    figure('position', [300 300 300 300])
    histVals =  [onMidgHist; offMidgHist; onParHist; offParHist];
    bar(histVector, histVals', 'grouped')
    set(gca, 'xlim', [0 0.5])
    xlabel('SD/mean of cumulative Gaussian fit')
    ylabel('number of cells')
    
    
    %overlapped (need to add transparencies)
    figure('position', [300 300 300 300]); hold on
    onParH = bar(histVector, onParHist');
    offParH = bar(histVector, offParHist');
    onMidgH = bar(histVector, onMidgHist');
    offMidgH = bar(histVector, offMidgHist');

    set(gca, 'xlim', [0 0.5])
    set(onParH,   'faceColor', colors(1,:), 'edgeColor', 'none')
    set(offParH,  'faceColor', colors(2,:), 'edgeColor', 'none')
    set(onMidgH,  'faceColor', colors(3,:), 'edgeColor', 'none')
    set(offMidgH, 'faceColor', colors(4,:), 'edgeColor', 'none')

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
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', colors(3,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'offMidg')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', colors(4,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', colors(1,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'offPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', colors(2,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigAllElecs, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
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
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', colors(3,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'offMidg')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', colors(4,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', colors(1,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'offPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', colors(2,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).maxSigStimElec, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
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
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', colors(3,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'offMidg')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', colors(4,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', colors(1,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        elseif strcmpi(cellInfo(ii).type, 'offPar')
            if cellInfo(ii).verMin
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', colors(2,:), 'MarkerSize', 6)
            else
                plot(cellInfo(ii).stimDistToCOM, cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
            end
        else
            error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
        end
    end
    hold off
    xlabel(['distance between stimulating electrode and', 10, 'electrical image center of mass (\mum)'], 'fontsize', 20)
    ylabel('threshold (µA)', 'fontsize', 20)
    set(gca, 'xlim', [-10 160], 'ylim', [0 5], 'xtick', [0 50 100 150], 'fontsize', 18)
    
    
    %%% plot thresh vs distance btwn stim elec and different offsets from COM in axon direction %%%
    
    %approximated angle
    axonAngle = 275; %in getElectrodeCoords61() coordinates: elec 13 on +x axis and elec 62 on +y axis
    
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
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', colors(3,:), 'MarkerSize', 6)
                else
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
                end
            elseif strcmpi(cellInfo(ii).type, 'offMidg')
                if cellInfo(ii).verMin
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', colors(4,:), 'MarkerSize', 6)
                else
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
                end
                
            elseif strcmpi(cellInfo(ii).type, 'onPar')
                if cellInfo(ii).verMin
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', colors(1,:), 'MarkerSize', 6)
                else
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
                end
            elseif strcmpi(cellInfo(ii).type, 'offPar')
                if cellInfo(ii).verMin
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', colors(2,:), 'MarkerSize', 6)
                else
                    plot(distsToOffsetCOM(ii,jj), cellInfo(ii).thresh, 'o', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', 'none', 'MarkerSize', 6)
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
