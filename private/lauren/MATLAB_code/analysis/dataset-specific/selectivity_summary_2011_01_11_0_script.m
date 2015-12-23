%%% note: all midget and parasol cells over array except purely axonal
%%% signals and cell 106 (essentially axonal signal, too small to
%%% analyze) are CHECKED for nontarget stimulation.

%%% However, only cells that are considered "over the array" (somatic
%%% signal meets criteria from removeSmallSigEdgeCells.m) are targeted
%%% for selective stimulation


% need to choose specific criteria for classification of nontarget
% stimulation as somatic or axonal

clear all

plotVsCurrent = true;

%% generates figure summarizing selectivity results for piece 2011-01-11-0

pathToAnalysis = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/';
pathToEi = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030-lh/data030-lh.ei';

cutoffMult = 0.5;

noStimCells = [286 391 187 381 395 444 624 755 879]; %cells that don't make it past 20% 
%or nontarget cell analysis can't be performed past 20% of target cell activation

%for loop only used to allow code folding
for ii = 1
    ind = 0;
    
    % some cells just starting to also be stimmed at end; 91 is only
    % analyzable through ~0.6 (m59)
    ind = ind+1;
    cellInfo(ind).targetCell     = 91;
    cellInfo(ind).otherCells     = [62 111 152];  %cell with lowest threshold, other than target cell (0 signifies no other cell responds within analyzed region)
    cellInfo(ind).otherStimAxon  = [false false false]; %false means somatic, true means axonal
    cellInfo(ind).maxMovAnalyzed = 63; %highest movie number that could be analyzed for other cell activity
    cellInfo(ind).stimElec       = 7;
    cellInfo(ind).type           = 'onMidg';

    %
    ind = ind+1; %no detected stim
    cellInfo(ind).targetCell     = 286;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'onMidg';

    
    %some overlap with somatic stim of 305
    ind = ind+1;
    cellInfo(ind).targetCell     = 302; %#ok<*SAGROW>
    cellInfo(ind).otherCells     = 305;
    cellInfo(ind).otherStimAxon  = false;
    cellInfo(ind).maxMovAnalyzed = 51;
    cellInfo(ind).stimElec       = 21;
    cellInfo(ind).type           = 'onMidg';

    %no detected stim
    ind = ind+1;
    cellInfo(ind).targetCell     = 391;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'onMidg';

    % thresh overlapping axon stim events but all above threshold of n424
    ind = ind+1;
    cellInfo(ind).targetCell     = 424;
    cellInfo(ind).otherCells     = [407 18 888];
    cellInfo(ind).otherStimAxon  = [true true true];
    cellInfo(ind).maxMovAnalyzed = 39;
    cellInfo(ind).stimElec       = 29;
    cellInfo(ind).type           = 'onMidg';
    
    % not selective - several cells/axons also respond
    ind = ind+1;
    cellInfo(ind).targetCell     = 556;
    cellInfo(ind).otherCells     = [767 578 647 754 879 573];
    cellInfo(ind).otherStimAxon  = [true false false true true false];
    cellInfo(ind).maxMovAnalyzed = 63;
    cellInfo(ind).stimElec       = 38;
    cellInfo(ind).type           = 'onMidg';

    % selective within tested current range (gets to 0.66)
    ind = ind+1;
    cellInfo(ind).targetCell     = 586;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 65;
    cellInfo(ind).stimElec       = 37;
    cellInfo(ind).type           = 'onMidg';
    
    % nonselective: axon stim and somatic stim at lower threshold than n631
    ind = ind+1;
    cellInfo(ind).targetCell     = 631;
    cellInfo(ind).otherCells     = [680 768 573];
    cellInfo(ind).otherStimAxon  = [false true false];
    cellInfo(ind).maxMovAnalyzed = 57;
    cellInfo(ind).stimElec       = 43;
    cellInfo(ind).type           = 'onMidg';
    
    % nonselective: 2 somatic stim events above n767 threshold
    ind = ind+1;
    cellInfo(ind).targetCell     = 767;
    cellInfo(ind).otherCells     = [754 768];
    cellInfo(ind).otherStimAxon  = [false false];
    cellInfo(ind).maxMovAnalyzed = 31;
    cellInfo(ind).stimElec       = 52;
    cellInfo(ind).type           = 'onMidg';
    
    % fully selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 856;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 47;
    cellInfo(ind).stimElec       = 15;
    cellInfo(ind).type           = 'onMidg';
    
    % fully selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 902;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 29;
    cellInfo(ind).stimElec       = 61;
    cellInfo(ind).type           = 'onMidg';
    
    % selective within analyzable region
    ind = ind+1;
    cellInfo(ind).targetCell     = 48;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 57;
    cellInfo(ind).stimElec       = 4;
    cellInfo(ind).type           = 'offMidg';
    
    % fully selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 62;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 65;
    cellInfo(ind).stimElec       = 5;
    cellInfo(ind).type           = 'offMidg';
    
    % fully selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 111;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 31;
    cellInfo(ind).stimElec       = 8;
    cellInfo(ind).type           = 'offMidg';
    
    % mostly selective: one axon stim event above n153 threshold
    ind = ind+1;
    cellInfo(ind).targetCell     = 153;
    cellInfo(ind).otherCells     = 62;
    cellInfo(ind).otherStimAxon  = true;
    cellInfo(ind).maxMovAnalyzed = 35;
    cellInfo(ind).stimElec       = 14;
    cellInfo(ind).type           = 'offMidg';
    
    % completely selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 169;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 33;
    cellInfo(ind).stimElec       = 12;
    cellInfo(ind).type           = 'offMidg';
    
    % no detectible stimulation
    ind = ind+1;
    cellInfo(ind).targetCell     = 187;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'offMidg';
    
    % almost completely selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 259;
    cellInfo(ind).otherCells     = 153;
    cellInfo(ind).otherStimAxon  = false;
    cellInfo(ind).maxMovAnalyzed = 33;
    cellInfo(ind).stimElec       = 18;
    cellInfo(ind).type           = 'offMidg';
    
    % nonselective:
    ind = ind+1;
    cellInfo(ind).targetCell     = 305;
    cellInfo(ind).otherCells     = [302 856];
    cellInfo(ind).otherStimAxon  = [false true];
    cellInfo(ind).maxMovAnalyzed = 61;
    cellInfo(ind).stimElec       = 21;
    cellInfo(ind).type           = 'offMidg';
    
    % one stim event below thresh of n338
    ind = ind+1;
    cellInfo(ind).targetCell     = 338;
    cellInfo(ind).otherCells     = 302;
    cellInfo(ind).otherStimAxon  = false;
    cellInfo(ind).maxMovAnalyzed = 45;
    cellInfo(ind).stimElec       = 23;
    cellInfo(ind).type           = 'offMidg';
    
    % not analyzable
    ind = ind+1;
    cellInfo(ind).targetCell     = 381;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'offMidg';
    
    % not analyzable
    ind = ind+1;
    cellInfo(ind).targetCell     = 395;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'offMidg';
    
    % no detectible stim
    ind = ind+1;
    cellInfo(ind).targetCell     = 444;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'offMidg';
    
    % fully selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 457;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 41;
    cellInfo(ind).stimElec       = 31;
    cellInfo(ind).type           = 'offMidg';
    
    % nonselective: multiple axon stimulation events and a somatic stim
    % event (n381 and n578 semi-axonal stim events)
    ind = ind+1;
    cellInfo(ind).targetCell     = 470;
    cellInfo(ind).otherCells     = [767 902 381 457 578 814];
    cellInfo(ind).otherStimAxon  = [true true false false false true];
    cellInfo(ind).maxMovAnalyzed = 55;
    cellInfo(ind).stimElec       = 31;
    cellInfo(ind).type           = 'offMidg';
    
    % nonselective: one somatic stim even at lower thresh
    ind = ind+1;
    cellInfo(ind).targetCell     = 578;
    cellInfo(ind).otherCells     = 573;
    cellInfo(ind).otherStimAxon  = false;
    cellInfo(ind).maxMovAnalyzed = 61;
    cellInfo(ind).stimElec       = 39;
    cellInfo(ind).type           = 'offMidg';
    
    % no detectible stim
    ind = ind+1;
    cellInfo(ind).targetCell     = 624;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'offMidg';
    
    % nonselective: one somatic stim and 3 axon stim below n647 threshold
    ind = ind+1;
    cellInfo(ind).targetCell     = 647;
    cellInfo(ind).otherCells     = [556 767 754 879];
    cellInfo(ind).otherStimAxon  = [false true true true];
    cellInfo(ind).maxMovAnalyzed = 63;
    cellInfo(ind).stimElec       = 38;
    cellInfo(ind).type           = 'offMidg';
    
    % nonselective: single axon stimulated
    ind = ind+1;
    cellInfo(ind).targetCell     = 680;
    cellInfo(ind).otherCells     = 768;
    cellInfo(ind).otherStimAxon  = true;
    cellInfo(ind).maxMovAnalyzed = 43;
    cellInfo(ind).stimElec       = 46;
    cellInfo(ind).type           = 'offMidg';
    
    % nonselective: 2 soma stimulated at lower than thresh
    ind = ind+1;
    cellInfo(ind).targetCell     = 754;
    cellInfo(ind).otherCells     = [767 768];
    cellInfo(ind).otherStimAxon  = [false false];
    cellInfo(ind).maxMovAnalyzed = 41;
    cellInfo(ind).stimElec       = 52;
    cellInfo(ind).type           = 'offMidg';
    
    % no detectible stim
    ind = ind+1;
    cellInfo(ind).targetCell     = 755;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'offMidg';
    
    % completely selective (figure cell)
    ind = ind+1;
    cellInfo(ind).targetCell     = 799;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 19;
    cellInfo(ind).stimElec       = 54;
    cellInfo(ind).type           = 'offMidg';
    
    % fully selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 814;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 47;
    cellInfo(ind).stimElec       = 49;
    cellInfo(ind).type           = 'offMidg';
    
    % fully selective: axon stim of n18 shouldn't actually overlap with
    % thresh +/- 2SD region of n877
    ind = ind+1;
    cellInfo(ind).targetCell     = 877;
    cellInfo(ind).otherCells     = 18;
    cellInfo(ind).otherStimAxon  = true;
    cellInfo(ind).maxMovAnalyzed = 37;
    cellInfo(ind).stimElec       = 47;
    cellInfo(ind).type           = 'offMidg';
    
    % no detectible stim
    ind = ind+1;
    cellInfo(ind).targetCell     = 879;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = [];
    cellInfo(ind).stimElec       = [];
    cellInfo(ind).type           = 'offMidg';
    
    % fully selective
    ind = ind+1;
    cellInfo(ind).targetCell     = 951;
    cellInfo(ind).otherCells     = [];
    cellInfo(ind).otherStimAxon  = [];
    cellInfo(ind).maxMovAnalyzed = 33;
    cellInfo(ind).stimElec       = 6;
    cellInfo(ind).type           = 'offMidg';
end


%remove cells that don't meet criteria to be "over the array"
[cellInfo cellInfoOffArray] = removeSmallSigEdgeCells(cellInfo, pathToEi, 'cutoffMult', cutoffMult);

%removes cells from 'noStimCells' list from struct array and counts number
%of cells that are "over the array" but still aren't stimulated or can't be
%fully analyzed
nCellsNoStimOverArray = 0;
for ii = length(cellInfo):-1:1
    if any(noStimCells == cellInfo(ii).targetCell)
        cellInfo(ii) = [];
        nCellsNoStimOverArray = nCellsNoStimOverArray+1;
    end
end

nCells = length(cellInfo);


%% get necessary info from elecResp files

normInv30 = norminv(0.8, 0, 1); %distance from 0.5 to 0.8 probability, in SDs
normInv45 = norminv(0.95, 0, 1);%distance from 0.5 to 0.95 probability, in SDs

for ii = 1:nCells
    load([pathToAnalysis 'elecResp_n' num2str(cellInfo(ii).targetCell) '_p' num2str(cellInfo(ii).stimElec) '_w50'])
    
    elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', 0); %bootstrapReps = 100
    
    cellInfo(ii).thresh = elecResp.analysis.threshold;
    cellInfo(ii).sigma = elecResp.analysis.threshold/(-sqrt(2)*elecResp.analysis.erfParams(2));
    
    %amplitudes at which response curve crosses 0.05, 0.2, 0.5, 0.8 and 0.95
    %cellInfo(ii).cdfAmps = cellInfo(ii).thresh + cellInfo(ii).sigma*[-normInv45 -normInv30 0 normInv30 normInv45];
    
    cellInfo(ii).maxAmpAnalyzed = abs(elecResp.stimInfo.stimAmps(elecResp.stimInfo.movieNos ==  cellInfo(ii).maxMovAnalyzed));
    
    %convert values to units of standard deviation relative to threshold
    cellInfo(ii).maxAnalyzedInSD = (cellInfo(ii).maxAmpAnalyzed - cellInfo(ii).thresh)/cellInfo(ii).sigma;
    
    %save([pathToAnalysis 'elecResp_n' num2str(cellInfo(ii).targetCell) '_p' num2str(cellInfo(ii).stimElec) '_w50'], 'elecResp')
    
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
    
    clear elecResp
    
    cellInfo(ii).otherThresh = zeros(size(cellInfo(ii).otherCells));
    cellInfo(ii).otherSigma = zeros(size(cellInfo(ii).otherCells));
    
    for jj = 1:length(cellInfo(ii).otherCells)
        load([pathToAnalysis 'elecResp_n' num2str(cellInfo(ii).otherCells(jj)) '_p' num2str(cellInfo(ii).stimElec) '_w50'])
        
        elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', 0); %bootstrapReps = 100
        
       
        %check that non-target cell was actually analyzed in the movie
        %that's supposedly the last movie analyzed, or else already passed
        %its threshold
        maxMovAnalyzedInd = find(elecResp.stimInfo.movieNos == cellInfo(ii).maxMovAnalyzed);
        maxMovAnalyzedLog = false;
        for kk = maxMovAnalyzedInd:length(elecResp.analysis.type)
            if ~isempty(elecResp.analysis.type{kk})
                maxMovAnalyzedLog = true;
                break
            end
        end
        if ~maxMovAnalyzedLog && max(elecResp.analysis.successRates) < 0.5
            error(['non-target cell ' num2str(cellInfo(ii).otherCells(jj)) ' was not analyzed in movie = maxMovAnalyzed',...
                10 '(target cell = ' num2str(cellInfo(ii).targetCell) ')'])
        end
        
        %make sure non-target cell passed threshold-SD and if it did, that
        %it also passed p = 0.2 (required for threshold estimate)
        if max(elecResp.analysis.successRates) < 0.155
            disp(['***non-target cell ' num2str(cellInfo(ii).otherCells(jj)) ' was not stimulated beyond 1 SD (p=0.155) from thresh***'...
                10 '(target cell = ' num2str(cellInfo(ii).targetCell) ')'])
        elseif max(elecResp.analysis.successRates) < 0.2 || min(elecResp.analysis.successRates > 0.8)
            error(['not enough of curve is measured for cell ' num2str(cellInfo(ii).otherCells(jj)) ' to estimate threshold'])
        end
        
        cellInfo(ii).otherThresh(jj) = elecResp.analysis.threshold;
        cellInfo(ii).otherSigma(jj) = 1/(sqrt(2)*elecResp.analysis.erfParams(1));
        
        %save([pathToAnalysis 'elecResp_n' num2str(cellInfo(ii).otherCells(jj)) '_p' num2str(cellInfo(ii).stimElec) '_w50'], 'elecResp')
        
        clear elecResp
    end
end

%% apply maxSuccRateCutOff; add cells that don't pass to "no stim" category
maxSuccRateCutoff = 0.4;

count = 0;
for ii = length(cellInfo):-1:1
    if max(cellInfo(ii).data(2,:)) < maxSuccRateCutoff
        disp(['cell ' num2str(cellInfo(ii).targetCell) ' didn''t meet success rate requirement of at least ' num2str(maxSuccRateCutoff)])
        
        %add cell to "cellInfoNoStim"
        count = count+1;
        cellInfoNoStim(count) = cellInfo(ii);
        
        %remove cell from "cellInfoNoStim"
        cellInfo(ii) = [];
    end
end

nCells = length(cellInfo);

%% determine max current amplitude with verified selectivity for each cell
for ii = 1:nCells
    maxSel = cellInfo(ii).maxAmpAnalyzed; %maximum amp in verified region
    for jj = 1:length(cellInfo(ii).otherCells)
        maxSel = min([maxSel, cellInfo(ii).otherThresh(jj) - normInv30*cellInfo(ii).otherSigma(jj)]);
    end
        
    cellInfo(ii).maxSelInSD = (maxSel - cellInfo(ii).thresh)/cellInfo(ii).sigma;
end

%% generate figure
figure('position', [100 0 400 900])
axes('position', [0.1 0.05 0.8 0.9])
hold on

% sort by size of analyzable region
%analyzedLims = [cellInfo.maxAnalyzedInSD];
%[sortedLims plotOrder] = sort(analyzedLims, 2, 'descend');

% in order to sort within cell type
onMidgBin = false(1, length(cellInfo));
for ii = 1:nCells
    if strcmpi(cellInfo(ii).type, 'onMidg')
        onMidgBin(ii) = true;
    end
end

if ~plotVsCurrent
    %order cells by type and selectivity (plotting goes from bottom to top)
    maxSelective = [cellInfo.maxSelInSD];
    [~, plotOrderOnMidg] = sort(maxSelective(onMidgBin), 2, 'descend');
    [~, plotOrderOffMidg] = sort(maxSelective(~onMidgBin), 2, 'descend');
    plotOrder = [plotOrderOffMidg+length(plotOrderOnMidg) 0 plotOrderOnMidg]; %zero so that there is a gap between ON and OFF
    
    %remove repeated fully selective cells
    plotOrder([1:8 20]) = [];
    nPlottedRows = length(plotOrder);
else
    %order cells by target cell threshold
    targThreshs = [cellInfo.thresh];
    [~, plotOrderOnMidg] = sort(targThreshs(onMidgBin), 2, 'descend');
    [~, plotOrderOffMidg] = sort(targThreshs(~onMidgBin), 2, 'descend');
    plotOrder = [plotOrderOffMidg+length(plotOrderOnMidg) 0 plotOrderOnMidg]; %zero so that there is a gap between ON and OFF
end

%determine plot order of non-target cells
for ii = 1:nCells
    otherBarMin = cellInfo(ii).otherThresh - cellInfo(ii).otherSigma*normInv30;
    if ~plotVsCurrent
        %to be coded
    else
        [otherBarMinOrdered, otherCellOrder] = sort(otherBarMin,2, 'descend');
        cellInfo(ii).otherPlotOrder = otherCellOrder(otherBarMinOrdered < cellInfo(ii).maxAmpAnalyzed);
    end
end



if ~plotVsCurrent
    %divide each target cell region up for bar plots
    barH = 0.8/7; %height of each non-target cell bar
    %yPlotRegion = [0.4-(barH)*length(cellInfo(ii).otherCells)/2 0.4+(barH)*length(cellInfo(ii).otherCells)/2]; %center the group of bars within region
    %yPlotRegion = [0 0.8]
    
    for ii = 1:nCells
        plotPos = find(plotOrder==ii);
        
        if ~isempty(plotPos)
            
            analyzedLim = cellInfo(ii).maxAnalyzedInSD;
            
            % plot type 3: bar from p = 0.2 to 0.8 with triangle denoting analysis limit
            
            %fill([-normInv30*[1 1] min([analyzedLim normInv30])*[1 1]],...
            %    plotPos+[1-0.02 1-barH+0.02 1-barH+0.02 1-0.02], [0.5 0.5 0.5], 'edgeColor', 'none')
            plot([-normInv30 min([analyzedLim normInv30])], plotPos+1-0.5*barH*[1 1], 'k-', 'linewidth', 2)
            
            if analyzedLim >= 0
                plot(0, plotPos+1-0.5*barH, 'o', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])
                %plot([0 0], plotPos+[1-0.02 1-barH+0.02], '-', 'lineWidth', 1, 'color', [0 0 0])
            end
            
            if analyzedLim < normInv45
                plot(analyzedLim*[1 1], plotPos+[1 1], 'v', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])
            end
            
            for jj = 1:length(cellInfo(ii).otherCells)
                
                %determine amplitudes at which nontarget cell crosses p = 0.2 and 0.8
                curveRange = cellInfo(ii).otherThresh(jj) + cellInfo(ii).otherSigma(jj)*normInv30*[-1 1];
                
                %curveRange = cellInfo(ii).otherThresh(jj) + cellInfo(ii).otherSigma(jj)*[-1 1]; %threshold +/- 1 SD
                %converts to units of target cell SDs from target cell thresh
                curveRange = (curveRange - cellInfo(ii).thresh)/cellInfo(ii).sigma;
                otherThresh = (cellInfo(ii).otherThresh(jj) - cellInfo(ii).thresh)/cellInfo(ii).sigma;
                
                if cellInfo(ii).otherStimAxon(jj)
                    %for now, don't distinguish between somatic + axonal stim
                    fillColor = [1 0.5 0.5];
                    threshColor = [1 0 0];
                else
                    %fillColor = [0.5 0.5 1];
                    %threshColor = [0 0 1];
                    fillColor = [1 0.5 0.5];
                    threshColor = [1 0 0];
                end
                
                %yPlotVals = plotPos + [yPlotRegion(2)-0.02+barH yPlotRegion(2)+0.02] - (jj+1)*barH;
                %yPlotVals = [plotPos+0.38+0.8/6 plotPos+0.42] - 0.8*jj/6;
                %yPlotVals = [yPlotVals yPlotVals(2) yPlotVals(1)];
                yPlotVals = plotPos + [-0.02 -barH+0.02 -barH+0.02 -0.02] + 1 - jj*barH;
                
                %plot non-target cells
                if curveRange(2) <= -normInv45 %completely off left side of plot
                    plot(-normInv45-0.1, mean(yPlotVals), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', threshColor)
                    %fill(-normInv45+[-0.2 -0.2 -0.1 -0.1], yPlotVals, fillColor,...
                    %    'edgecolor', 'none')
                elseif curveRange(1) <= -normInv45 %partially off left side of plot
                    %fill([-normInv45-0.2 -normInv45-0.2 curveRange(2) curveRange(2)], yPlotVals, fillColor,...
                    %    'edgecolor', 'none')
                    if otherThresh < -normInv45
                        plot(-normInv45-0.1, mean(yPlotVals), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', threshColor)
                    end
                    %fill([-normInv45-0.1 -normInv45-0.1 curveRange(2) curveRange(2)], yPlotVals, fillColor,...
                    %    'edgecolor', 'none')
                    plot([-normInv45-0.1 curveRange(2)], mean(yPlotVals)*[1 1], 'r-', 'linewidth', 2)
                elseif curveRange(2) <= analyzedLim %fully on plot
                    %fill([curveRange(1) curveRange(1) curveRange(2) curveRange(2)], yPlotVals, fillColor,...
                    %    'edgecolor', 'none')
                    plot([curveRange(1) curveRange(2)], mean(yPlotVals)*[1 1], 'r-', 'linewidth', 2)
                    
                elseif curveRange(1) <= min([analyzedLim normInv45]) %partially off right side of plot
                    %fill([curveRange(1) curveRange(1) analyzedLim analyzedLim], yPlotVals, fillColor,...
                    %    'edgecolor', 'none')
                    plot([curveRange(1) analyzedLim], mean(yPlotVals)*[1 1], 'r-', 'linewidth', 2)
                    
                end
                if otherThresh <= analyzedLim && otherThresh >= -normInv45 %threshold within plot
                    plot(otherThresh, mean(yPlotVals), 'o', 'MarkerFaceColor', threshColor, 'MarkerEdgeColor', 'none')
                    %plot([otherThresh otherThresh], yPlotVals(1:2), '-', 'lineWidth', 2, 'color', threshColor)
                    %plot([otherThresh otherThresh], yPlotVals(1:2) + [0.2/6, -0.2/6], '-', 'lineWidth', 2, 'color', fillColor)
                end
            end
        end
    end
else %plot in terms of current
    plotPosCurrent = 0;
    for ii = 1:length(plotOrder)
        cellInd = plotOrder(ii);
        if cellInd ~= 0
            lim = cellInfo(cellInd).maxAmpAnalyzed;
            nCellsInRegion = length(cellInfo(cellInd).otherPlotOrder)+1;
            plotRegionX = [0 0 lim lim];
            plotRegionY = plotPosCurrent + 0.5 + [0 nCellsInRegion nCellsInRegion 0];
            fill(plotRegionX, plotRegionY, [0.8 0.8 0.8], 'edgecolor', 'none');
            
            otherCellLowerLim = 10;
            for jj = 1:length(cellInfo(cellInd).otherPlotOrder);
                plotPosCurrent = plotPosCurrent+1;
                sig = cellInfo(cellInd).otherSigma(cellInfo(cellInd).otherPlotOrder(jj));
                t = cellInfo(cellInd).otherThresh(cellInfo(cellInd).otherPlotOrder(jj));
                
                otherCellLowerLim = min([otherCellLowerLim, t-normInv30*sig]);
                                
                if t-normInv30*sig <= lim
                    plot([t-normInv30*sig min([lim t+normInv30*sig])], plotPosCurrent*[1 1], '-', 'linewidth', 2, 'color', [1 0 0])
                end
                
                if lim >= t
                    plot(t, plotPosCurrent, 'o', 'markerEdgeColor', 'none', 'markerFaceColor', [1 0 0])
                end
            end
            
            sig = cellInfo(cellInd).sigma;
            t = cellInfo(cellInd).thresh;
            
            plotPosCurrent = plotPosCurrent+1;
            plot([t-normInv30*sig min([lim t+normInv30*sig])], plotPosCurrent*[1 1], 'k-', 'linewidth', 2)
            if lim >= t
                plot(t, plotPosCurrent, 'o', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])
            end
            
            %mark which target cells are "selective"
%             if otherCellLowerLim >= t+normInv30*sig && lim >= t+normInv30*sig %verified non-overlapping
%                 plot(0.1, plotPosCurrent, 'k*')
%             elseif otherCellLowerLim >= t+normInv30*sig %unverified non-overlapping
%                 plot(0.1, plotPosCurrent, 'k+')
%             end
            
            %plot(lim*[1 1], plotPosCurrent+[1 1], 'v', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])
            
            
        else %put a little white space between cell types
            plotPosCurrent = plotPosCurrent+2;
        end
        plotPosCurrent = plotPosCurrent+1;
    end
end

%plot response curve to illustrate boundaries
%fill([-2 -2 2 2], nCells+[3 5 5 3], [0.9 0.9 0.9], 'edgeColor', 'none')
%fill([-1 -1 1 1], nCells+[3 5 5 3], [0.8 0.8 0.8], 'edgeColor', 'none')
%plot([0 0], nPlottedRows+[3 5], '-', 'lineWidth', 2, 'color', [0 0 0])

%xgrid = -normInv45:0.01:normInv45+0.01;
%erfCurve = 0.5 + 0.5*erf((1/sqrt(2))*xgrid(1,:));

%plot(xgrid, erfCurve*2+nPlottedRows+3, 'k-', 'linewidth', 2)

hold off
if ~plotVsCurrent
    plot([-normInv45 -normInv45], [0 nPlottedRows+2.5], 'k--')
    
    set(gca, 'xlim', [-normInv45-0.2 normInv45], 'ylim', [1 nPlottedRows+1.5], 'ytick', [], 'xtick', [-normInv45 -normInv30 0 normInv30 normInv45], 'xticklabel', [0.05 0.2 0.5 0.8 0.95])
    xlabel('target cell response probability')
else
    xlabel('current amplitude (µA)')
    set(gcf, 'color', [1 1 1], 'position', [100 0 400 8*plotPosCurrent+100])
    set(gca, 'ytick', [], 'ycolor', [1 1 1], 'ylim', [0 plotPosCurrent],...
        'xlim', [0 4.5], 'xtick', [0 1.5 3 4.5], 'units', 'pixels', 'position', [40 50 320 8*plotPosCurrent])
    
    %plot inset
    axes('units', 'pixels', 'position', [260 8*plotPosCurrent-30 100 80], 'box', 'on'); hold on
    plot(-2:0.002:2, normcdf(-2:0.002:2, 0, 1), '-', 'lineWidth', 2, 'color', [0.8 0.8 0.8])
    plot(-normInv30:0.002:normInv30, normcdf(-normInv30:0.002:normInv30, 0, 1), 'k-', 'lineWidth', 2)
    plot(0, 0.5, 'o', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])

    plot([-2 -normInv30], [0.2 0.2], 'k--')
    plot([-2 normInv30], [0.8 0.8], 'k--')
    plot(-normInv30*[1 1], [0.05 0.2], 'k--')
    plot(normInv30*[1 1], [0.05 0.8], 'k--')

    
    plot([-normInv30 normInv30], [0.05 0.05], 'k-', 'lineWidth', 2)
    plot(0, 0.05, 'o', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])

    set(gca, 'ytick', [0 0.2 0.5 0.8 1], 'xtick', [-normInv30 0 normInv30], 'xticklabel', [])
end

%% plot of selectivity vs. threshold

PW = 50;
figure('position', [100 100 800 400])

%probability of target cell
axes('position', [0.08 0.2 0.4 0.7]); hold on
for ii = 1:nCells
    if cellInfo(ii).maxAnalyzedInSD <= cellInfo(ii).maxSelInSD &&...
            cellInfo(ii).maxAnalyzedInSD < normInv45 %no other stim present but not full analysis
        plot(cellInfo(ii).thresh*PW*[1 1], normcdf(cellInfo(ii).maxSelInSD, 0, 1), 'k^')
    else %full analysis or other stim before analysis limit
        plot(cellInfo(ii).thresh*PW, normcdf(min([cellInfo(ii).maxSelInSD normInv45]), 0, 1), 'ko')

    end
end
xlabel('threshold charge amplitude (pA)')
ylabel('selectivity index (target cell response prob.)')


%amplitude, in target cell SDs from threshold
axes('position', [0.55 0.2 0.4 0.7]); hold on
for ii = 1:nCells
    if cellInfo(ii).maxAnalyzedInSD <= cellInfo(ii).maxSelInSD &&...
            cellInfo(ii).maxAnalyzedInSD < normInv45 %no other stim present but not full analysis
        plot(cellInfo(ii).thresh*PW*[1 1], cellInfo(ii).maxSelInSD, 'r^')
    else %full analysis or other stim before analysis limit
        plot(cellInfo(ii).thresh*PW, min([cellInfo(ii).maxSelInSD normInv45]), 'ro')
    end
end

xlabel('threshold charge amplitude (pA)')
ylabel('selectivity index (SDs from target cell thresh)')

%return

%% examples
keyboard
datarun = load_data('2011-01-11-0/data030-lh/data030-lh');
datarun = load_index(datarun);
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));
datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
%datarun = get_sta_fits_from_vision(datarun,'all');
datarun.names.obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030-lh';
datarun = load_obvius_sta_fits(datarun);
datarun = get_sta_fits_from_obvius(datarun, 'all');


%calculate electrode positions
datarun = compute_monitor_to_array_transformation(datarun);
T_array_to_monitor = fliptform(datarun.piece.T_monitor_to_array);
T_monitor_to_STA = fliptform(coordinate_transform(datarun, 'monitor'));
T = maketform('composite',[T_monitor_to_STA T_array_to_monitor]);
e_pos = tformfwd(T, datarun.ei.position);


%these lists should contain all cells were checked for nontarget activation, 
%which are all of the parasol and midget cells in vision analysis
%except those under "axon" or "off array" categories
goodCells{1} = [2 152 257 407 436 453 661 768 871]; %ON parasol
goodCells{2} = [18 138 288 392 573 617 751 753 888]; %OFF parasol
goodCells{3} = [91 286 302 391 424 556 586 631 767 856 902]; %ON midget
goodCells{4} = [48 62 111 153 169 187 259 305 338 381 395 444 457 470 578 624 647 680 754 755 799 814 877 879 951]; %OFF midget


axesWidthRatio = 2*normInv45/(2*normInv45+0.2); %so that axis lengths are equivalent with summary plot


%% 

plotColors = jet(101);
axesPositions = {[0 0 0.5 0.5], [0.5 0 0.5 0.5], [0 0.5 0.5 0.5], [0.5 0.5 0.5 0.5]};

fillColor = [0.5 0.5 1];
threshColor = [0 0 1];

%% cell 573, p39, m37 (fully selective, OFF parasol)

cellID = 573;
stimElec = 39;
movieNo = 37;
elecRespSuffix = '_w50';

selectivity_mosaic_plot(datarun, goodCells, stimElec, movieNo, pathToAnalysis,...
    'elecRespSuffix', elecRespSuffix)

%% cell 799, p54, m17 (fully selective)

cellID = 799;
stimElec = 54;
movie = 17;
movieNo = 17;
elecRespSuffix = '_w50';

selectivity_mosaic_plot(datarun, goodCells, stimElec, movieNo, pathToAnalysis,...
    'elecRespSuffix', elecRespSuffix)

%response curves plot
if true
    figure('position', [100 100 400 300])

    % response curves
    axes('position', [0.1 0.2 0.8*axesWidthRatio 0.7])

    hold on
    
    cellInd = find([cellInfo.targetCell]==cellID);
    plotXLim = cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*[-normInv45 normInv45];
    xTicks = cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*[-normInv45 -normInv30 0 normInv30 normInv45];

    %fill([plotXLim(1) plotXLim(1) plotXLim(2) plotXLim(2)], [0 1 1 0], [0.9 0.9 0.9], 'edgeColor', 'none')
    %fill(cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*[-1 -1 1 1], [0 1 1 0], [0.8 0.8 0.8], 'edgeColor', 'none')
    %plot(cellInfo(cellInd).thresh*[1 1], [0 1], 'color', [0.5 0.5 0.5], 'lineWidth', 2)
    
    xGrid = linspace(plotXLim(1), plotXLim(2), 1000);
    
    %plot target cell data
    load([pathToAnalysis 'elecResp_n' num2str(cellID) '_p' num2str(stimElec) '_w50'])
    
    data = zeros(3, length(elecResp.stimInfo.stimAmps));
    data(1,:) = abs(elecResp.stimInfo.stimAmps);
    data(2,:) = elecResp.analysis.successRates;
    for jj = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{jj})
            data(:,jj) = [];
        end
    end
    erfParams = elecResp.analysis.erfParams;
    mosaicStimAmp = abs(elecResp.stimInfo.stimAmps(elecResp.stimInfo.movieNos == movie));
    clear elecResp

    erfCurve = 0.5 + 0.5*erf(erfParams(1)*xGrid+erfParams(2));
    xGridInner = xGrid >= cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*-normInv30 &...
        xGrid <= cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*normInv30;
    
    %plot(xGrid, erfCurve, '-', 'lineWidth', 1, 'color', [0 0 0])
    plot(xGrid, erfCurve, '-', 'lineWidth', 2, 'color', [0.8 0.8 0.8])
    plot(xGrid(xGridInner), erfCurve(xGridInner), 'k-', 'lineWidth', 2)
    plot(data(1,:), data(2,:), 'o', 'markerEdgeColor', [0 0 0], 'markerFaceColor', [0 0 0])
    plot(data(1,data(1,:)==mosaicStimAmp), data(2,data(1,:)==mosaicStimAmp),...
            'o', 'markerEdgeColor', [0 0 0], 'markerFaceColor', [1 1 1])
        
    %plot nontarget cell data - none for this example
    
    plot(mosaicStimAmp*[1 1], [0 1], 'k--')
    
    set(gca, 'xlim', plotXLim, 'ylim', [0 1], 'ytick', [0 0.5 1])
    ylabel('response probability')
    xlabel('current amplitude (µA)')
    %xlabel('target cell response probability')

end


%% cell 302, p21, m47 (partially selective)

cellID = 302;
stimElec = 21;
movie = 47;
movieNo = 47;
elecRespSuffix = '_w50';


%mosaics plot

selectivity_mosaic_plot(datarun, goodCells, stimElec, movieNo, pathToAnalysis,...
    'elecRespSuffix', elecRespSuffix)


%response curves plot
if true
    figure('position', [100 100 400 300])
    axes('position', [0.1 0.2 0.8*axesWidthRatio 0.7])
    hold on
    
    cellInd = find([cellInfo.targetCell]==cellID);
    plotXLim = cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*[-normInv45 normInv45];

    
    %fill([plotXLim(1) plotXLim(1) plotXLim(2) plotXLim(2)], [0 1 1 0], [0.9 0.9 0.9], 'edgeColor', 'none')
    %fill(cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*[-1 -1 1 1], [0 1 1 0], [0.8 0.8 0.8], 'edgeColor', 'none')
    %plot(cellInfo(cellInd).thresh*[1 1], [0 1], 'color', [0.5 0.5 0.5], 'lineWidth', 2)
    
    xGrid = linspace(plotXLim(1), plotXLim(2), 1000);
    
    %plot target cell data
    load([pathToAnalysis 'elecResp_n' num2str(cellID) '_p' num2str(stimElec) '_w50'])
    
    data = zeros(3, length(elecResp.stimInfo.stimAmps));
    data(1,:) = abs(elecResp.stimInfo.stimAmps);
    data(2,:) = elecResp.analysis.successRates;
    for jj = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{jj})
            data(:,jj) = [];
        end
    end
    erfParams = elecResp.analysis.erfParams;
    mosaicStimAmp = abs(elecResp.stimInfo.stimAmps(elecResp.stimInfo.movieNos == movie));
    clear elecResp

    erfCurve = 0.5 + 0.5*erf(erfParams(1)*xGrid+erfParams(2));
    xGridInner = xGrid >= cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*-normInv30 &...
        xGrid <= cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*normInv30;
    
    plot(xGrid, erfCurve, '-', 'lineWidth', 2, 'color', [0.8 0.8 0.8])
    plot(xGrid(xGridInner), erfCurve(xGridInner), 'k-', 'lineWidth', 2)
    
    plot(data(1,:), data(2,:), 'o', 'markerEdgeColor', [0 0 0], 'markerFaceColor', [0 0 0])
    plot(data(1,data(1,:)==mosaicStimAmp), data(2,data(1,:)==mosaicStimAmp),...
            'o', 'markerEdgeColor', [0 0 0], 'markerFaceColor', [1 1 1])

    %plot nontarget cell data
    otherCells = cellInfo([cellInfo.targetCell]==cellID).otherCells;
    for kk = 1:length(otherCells)
        load([pathToAnalysis 'elecResp_n' num2str(otherCells(kk)) '_p' num2str(stimElec) '_w50'])
        data = zeros(3, length(elecResp.stimInfo.stimAmps));
        data(1,:) = abs(elecResp.stimInfo.stimAmps);
        data(2,:) = elecResp.analysis.successRates;
        for jj = length(elecResp.stimInfo.stimAmps): -1: 1
            if isempty(elecResp.analysis.type{jj})
                data(:,jj) = [];
            end
        end
        erfParams = elecResp.analysis.erfParams;
        clear elecResp
        
        erfCurve = 0.5 + 0.5*erf(erfParams(1)*xGrid+erfParams(2));
        xGridInner = xGrid >= cellInfo([cellInfo.targetCell]==cellID).otherThresh(kk) + cellInfo([cellInfo.targetCell]==cellID).otherSigma(kk)*-normInv30 &...
            xGrid <= cellInfo([cellInfo.targetCell]==cellID).otherThresh(kk) + cellInfo([cellInfo.targetCell]==cellID).otherSigma(kk)*normInv30;
        
        plot(xGrid, erfCurve, 'lineWidth', 2, 'color', [1 0.5 0.5])
        plot(xGrid(xGridInner), erfCurve(xGridInner), 'lineWidth', 2, 'color', [1 0 0])
        plot(data(1,:), data(2,:), 'o', 'markerEdgeColor', [1 0 0], 'markerFaceColor', [1 0 0])
        plot(data(1,data(1,:)==mosaicStimAmp), data(2,data(1,:)==mosaicStimAmp),...
            'o', 'markerEdgeColor', [1 0 0], 'markerFaceColor', [1 1 1])
    end
    plot(mosaicStimAmp*[1 1], [0 1], 'k--')

    set(gca, 'xlim', plotXLim, 'ylim', [0 1], 'ytick', [0 0.5 1])
    ylabel('response probability')
    xlabel('current amplitude (µA)')
    
%     %summary
%     axes('position', [0.1 0.1 0.8 0.05])
%     hold on
%     fill([plotXLim(1) plotXLim(1) plotXLim(2) plotXLim(2)], [0 1 1 0], [0.9 0.9 0.9], 'edgeColor', 'none')
%     fill(cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*[-1 -1 1 1], [0 1 1 0], [0.8 0.8 0.8], 'edgeColor', 'none')
%     plot(cellInfo(cellInd).thresh*[1 1], [0 1], 'color', [0.5 0.5 0.5], 'lineWidth', 2)
%     
%     for kk = 1:length(otherCells)
%         curveRange = [cellInfo(cellInd).otherThresh(kk)-cellInfo(cellInd).otherSigma(kk)...
%             cellInfo(cellInd).otherThresh(kk)+cellInfo(cellInd).otherSigma(kk)]; %threshold +/- 1 SD
%         fill([curveRange(1) curveRange(1) curveRange(2) curveRange(2)], [0.375 0.625 0.625 0.375], fillColor,...
%             'edgecolor', 'none')
%         plot(cellInfo(cellInd).otherThresh(kk)*[1 1], [0.375 0.625], '-', 'lineWidth', 2, 'color', threshColor)
%     end
%     set(gca, 'xlim', plotXLim, 'ylim', [0 1], 'ytick', [])
%     xlabel('current amplitude (µA)')
end


%% cell 631, p43, m51 (unselective; can't directly determine
% response rate of n768 at this amplitude, so extimated from curve fit)

cellID = 631;
stimElec = 43;
movie = 51;
movieNo = 51;
elecRespSuffix = '_w50';

%estimate success rate for n768 from response curve fit
load([pathToAnalysis 'elecResp_n768_p' num2str(stimElec) elecRespSuffix])
movieInd = find(elecResp.stimInfo.movieNos == movieNo);
succRate = 0.5 + 0.5*erf(elecResp.analysis.erfParams(1)*abs(elecResp.stimInfo.stimAmps(movieInd))+elecResp.analysis.erfParams(2));
succRateOverride = [768 succRate];

selectivity_mosaic_plot(datarun, goodCells, stimElec, movieNo, pathToAnalysis,...
    'elecRespSuffix', elecRespSuffix, 'succRateOverride', succRateOverride)

%response curves plot
if true
    figure('position', [100 100 400 300])
    axes('position', [0.1 0.2 0.8*axesWidthRatio 0.7])
    hold on
    
    cellInd = find([cellInfo.targetCell]==cellID);
    plotXLim = cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*[-normInv45 normInv45];
    
    maxAnalyzed = cellInfo(cellInd).maxAmpAnalyzed;
    
    %fill([plotXLim(1) plotXLim(1) maxAnalyzed maxAnalyzed], [0 1 1 0], [0.9 0.9 0.9], 'edgeColor', 'none')
    %fill(cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*[-1 -1 1 1], [0 1 1 0], [0.8 0.8 0.8], 'edgeColor', 'none')
    %plot(cellInfo(cellInd).thresh*[1 1], [0 1], 'color', [0.5 0.5 0.5], 'lineWidth', 2)
    plot(maxAnalyzed, 1, 'v', 'markerFaceColor', [0 0 0], 'markerEdgeColor', 'none')
    
    xGrid = linspace(plotXLim(1), maxAnalyzed, 1000);
    
    %plot target cell data
    load([pathToAnalysis 'elecResp_n' num2str(cellID) '_p' num2str(stimElec) '_w50'])
    
    data = zeros(3, length(elecResp.stimInfo.stimAmps));
    data(1,:) = abs(elecResp.stimInfo.stimAmps);
    data(2,:) = elecResp.analysis.successRates;
    for jj = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{jj})
            data(:,jj) = [];
        end
    end
    erfParams = elecResp.analysis.erfParams;
    mosaicStimAmp = abs(elecResp.stimInfo.stimAmps(elecResp.stimInfo.movieNos == movie));
    clear elecResp
    
    
    erfCurve = 0.5 + 0.5*erf(erfParams(1)*xGrid+erfParams(2));
    xGridInner = xGrid >= cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*-normInv30 &...
        xGrid <= cellInfo(cellInd).thresh + cellInfo(cellInd).sigma*normInv30;
    
    plot(xGrid, erfCurve, '-', 'lineWidth', 2, 'color', [0.8 0.8 0.8])
    plot(xGrid(xGridInner), erfCurve(xGridInner), 'k-', 'lineWidth', 2)
    
    plot(data(1,:), data(2,:), 'o', 'markerEdgeColor', [0 0 0], 'markerFaceColor', [0 0 0])
    plot(data(1,data(1,:)==mosaicStimAmp), data(2,data(1,:)==mosaicStimAmp),...
        'o', 'markerEdgeColor', [0 0 0], 'markerFaceColor', [1 1 1])

    %plot nontarget cell data
    otherCells = cellInfo([cellInfo.targetCell]==cellID).otherCells;
    for kk = 1:length(otherCells)
        load([pathToAnalysis 'elecResp_n' num2str(otherCells(kk)) '_p' num2str(stimElec) '_w50'])
        data = zeros(3, length(elecResp.stimInfo.stimAmps));
        data(1,:) = abs(elecResp.stimInfo.stimAmps);
        data(2,:) = elecResp.analysis.successRates;
        for jj = length(elecResp.stimInfo.stimAmps): -1: 1
            if isempty(elecResp.analysis.type{jj})
                data(:,jj) = [];
            end
        end
        erfParams = elecResp.analysis.erfParams;
        clear elecResp
        
        erfCurve = 0.5 + 0.5*erf(erfParams(1)*xGrid+erfParams(2));
        xGridInner = xGrid >= cellInfo([cellInfo.targetCell]==cellID).otherThresh(kk) + cellInfo([cellInfo.targetCell]==cellID).otherSigma(kk)*-normInv30 &...
            xGrid <= cellInfo([cellInfo.targetCell]==cellID).otherThresh(kk) + cellInfo([cellInfo.targetCell]==cellID).otherSigma(kk)*normInv30;
        
        plot(xGrid, erfCurve, 'lineWidth', 2, 'color', [1 0.5 0.5])
        plot(xGrid(xGridInner), erfCurve(xGridInner), 'lineWidth', 2, 'color', [1 0 0])
        
        plot(data(1,data(1,:)<=maxAnalyzed), data(2,data(1,:)<=maxAnalyzed), 'o', 'markerEdgeColor', [1 0 0], 'markerFaceColor', [1 0 0])
        plot(data(1,data(1,:)==mosaicStimAmp), data(2,data(1,:)==mosaicStimAmp),...
            'o', 'markerEdgeColor', [1 0 0], 'markerFaceColor', [1 1 1])
    end
    plot(mosaicStimAmp*[1 1], [0 1], 'k--')

    set(gca, 'xlim', plotXLim, 'ylim', [0 1], 'ytick', [0 0.5 1])
    ylabel('response probability')
    xlabel('current amplitude (µA)')
    
end


%% generates colorbar for figure

nColors = size(plotColors,1);

%plot horizontal color scale bar
figure('position', [100 100 500 100])
axes('position', [0.1 0.1 0.8 0.4])
hold on

yPoly = [0 1 1 0];
for i = 1:nColors
    xPoly = [(i-1)/nColors (i-1)/nColors i/nColors i/nColors];
    fill(xPoly, yPoly, plotColors(i,:), 'EdgeColor', 'none')
end

text(1, 1.5, '1', 'HorizontalAlignment', 'center', 'fontsize', 14)
text(1/nColors, 1.5, '0', 'HorizontalAlignment', 'center', 'fontsize', 14)
text(0.5*(nColors+1)/nColors, 1.5, '0.5', 'HorizontalAlignment', 'center', 'fontsize', 14)

% for i = 1:length(colorBarLabels)
%     if colorBarLabels(i) > zMin && colorBarLabels(i) < zMax
%         y = (colorBarLabels(i)-zMin)/(zMax-zMin);
%         plot([1 1.1], [y y], 'k-')
%         text(1.2, y, [num2str(colorBarLabels(i)), '%'])
%     end
% end
set(gca, 'ylim', [0 2], 'xlim', [0 1])
axis off
hold off


%% tally non-target activations that are putative axonal stimulation (stim
% elec > 150 µm from peak EI signal elec)

nNonTargetSomatic = 0;
nNonTargetAxonal = 0;
axonDistThresh = 150; %distance between stim elec and max EI signal that must be exceeded to consider it "axonal stimulation"
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
[xCoords yCoords] = getElectrodeCoords61();

for ii = 1:nCells
    analyzedLim = cellInfo(ii).maxAnalyzedInSD;
    
    %text(2, plotPos, ['n' num2str(cellInfo(ii).targetCell) ' p' num2str(cellInfo(ii).stimElec)])
    
    for jj = 1:length(cellInfo(ii).otherCells)
        
        curveRange = cellInfo(ii).otherThresh(jj) + cellInfo(ii).otherSigma(jj)*normInv30*[-1 1];
        
        %converts to units of target cell SDs from target cell thresh
        curveRange = (curveRange - cellInfo(ii).thresh)/cellInfo(ii).sigma;
        
        if curveRange(1) < min([normInv45 cellInfo(ii).maxAnalyzedInSD]) %non-target stim appears in selectivity summary plot
            
            %determine electrode with max EI signal
            ei = eiFile.getImage(cellInfo(ii).otherCells(jj));
            %calculate negative peak (absolute) signal on each electrode
            ei = squeeze(ei(1,2:end,:));
            eiAmps = max(-ei, [], 2);
            
            maxSigRecElec = find(eiAmps == max(eiAmps));
            
            stimDist = 30*norm([xCoords(maxSigRecElec) - xCoords(cellInfo(ii).stimElec) yCoords(maxSigRecElec) - yCoords(cellInfo(ii).stimElec)]);
            
            if stimDist > axonDistThresh
                nNonTargetAxonal =  nNonTargetAxonal + 1;
                if ~cellInfo(ii).otherStimAxon(jj) %make sure calculations are correct
                    keyboard
                end
            else
                nNonTargetSomatic = nNonTargetSomatic + 1;
                if cellInfo(ii).otherStimAxon(jj) %make sure calculations are correct
                    keyboard
                end
            end
        end
        
    end
end
