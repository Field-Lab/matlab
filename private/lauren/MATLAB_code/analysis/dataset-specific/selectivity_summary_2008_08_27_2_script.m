%% generates figure summarizing selectivity results for piece 2008-08-27-2
% checked for accuracy 2012-01-09

%%% note: all midget and parasol cells over array except one cell that 
%%% is too small to analyze (cell 437) are CHECKED for nontarget stimulation.

%%% However, only cells that are considered "over the array" (somatic
%%% signal meets criteria from removeSmallSigEdgeCells.m) are targeted
%%% for selective stimulation

clear all

plotVsCurrent = true;

pathToAnalysis = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/';
pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';

cutoffMult = 0.5;

noStimCells = [138 166 241 437 496 526 527 590 707 857]; %cells that don't make it past 20%

%for loop only used to allow code folding
for ii = 1
    % only gets to 0.4
    cellInfo(1).targetCell     = 18;
    cellInfo(1).otherCells     = [];  %cell with lowest threshold, other than target cell (0 signifies no other cell responds within analyzed region)
    cellInfo(1).otherStimAxon  = []; %false means somatic, true means axonal
    cellInfo(1).maxMovAnalyzed = 26; %highest movie number that could be analyzed for other cell activity
    cellInfo(1).stimElec       = 2;
    cellInfo(1).type           = 'onMidg';

    
    % selective with full curve but on edge of mosaic
    cellInfo(2).targetCell     = 91;
    cellInfo(2).otherCells     = 92;
    cellInfo(2).otherStimAxon  = false;
    cellInfo(2).maxMovAnalyzed = 26;
    cellInfo(2).stimElec       = 7;
    cellInfo(2).type           = 'onMidg';


    % elec 6: local min thresh, n931 has slightly lower thresh
    % previously elec 3: only gets to 0.5, and thresh (from fitting) is > than final amplitude
    cellInfo(3).targetCell     = 92;
    cellInfo(3).otherCells     = 931;
    cellInfo(3).otherStimAxon  = true;
    cellInfo(3).maxMovAnalyzed = 26;
    cellInfo(3).stimElec       = 6;
    cellInfo(3).type           = 'onMidg';

    %can't analyze completely (past m23) and no sign of stimulation, corner of
    %mosaic and edge of array
    cellInfo(4).targetCell     = 138;
    cellInfo(4).otherCells     = [];
    cellInfo(4).otherStimAxon  = [];
    cellInfo(4).maxMovAnalyzed = 23;
    cellInfo(4).stimElec       = [];
    cellInfo(4).type           = 'onMidg';

    %no sign of stimulation, but neighbors electrode 4, and can't analyze past
    %m25
    cellInfo(5).targetCell     = 166;
    cellInfo(5).otherCells     = [];
    cellInfo(5).otherStimAxon  = [];
    cellInfo(5).maxMovAnalyzed = 25;
    cellInfo(5).stimElec       = [];
    cellInfo(5).type           = 'onMidg';

    %completely selective; can analyze further if needed
    cellInfo(6).targetCell     = 167;
    cellInfo(6).otherCells     = [];
    cellInfo(6).otherStimAxon  = [];
    cellInfo(6).maxMovAnalyzed = 20;
    cellInfo(6).stimElec       = 12;
    cellInfo(6).type           = 'onMidg';

    %selective but only gets to 0.5 within measured range
    cellInfo(7).targetCell     = 182;
    cellInfo(7).otherCells     = [];
    cellInfo(7).otherStimAxon  = [];
    cellInfo(7).maxMovAnalyzed = 26;
    cellInfo(7).stimElec       = 13;
    cellInfo(7).type           = 'onMidg';

    %no sign of stimulation above 0.16, can only measure up to m24 on some
    %other electrodes (no sign of stim there either)
    cellInfo(8).targetCell     = 241;
    cellInfo(8).otherCells     = [];
    cellInfo(8).otherStimAxon  = [];
    cellInfo(8).maxMovAnalyzed = [];
    cellInfo(8).stimElec       = [];
    cellInfo(8).type           = 'onMidg';

    %unselective (axon) and only gets to 0.3; can't analyze 2 of the neighboring possible stim electrodes
    cellInfo(9).targetCell     = 256;
    cellInfo(9).otherCells     = [92 886];
    cellInfo(9).otherStimAxon  = [true true];
    cellInfo(9).maxMovAnalyzed = 26;
    cellInfo(9).stimElec       = 20;
    cellInfo(9).type           = 'onMidg';

    %nonselective, but stim electrode on corner of array
    cellInfo(10).targetCell     = 274;
    cellInfo(10).otherCells     = [257 902];
    cellInfo(10).otherStimAxon  = [false true];
    cellInfo(10).maxMovAnalyzed = 26;
    cellInfo(10).stimElec       = 19;
    cellInfo(10).type           = 'onMidg';

    %axonal stim above threshold
    cellInfo(11).targetCell     = 332;
    cellInfo(11).otherCells     = [782 811];
    cellInfo(11).otherStimAxon  = [true true];
    %cellInfo(11).otherCells     = [317 782 811]; %317 doesn't show up on plot (above analyzed region)
    %cellInfo(11).otherStimAxon  = [true true true];
    cellInfo(11).maxMovAnalyzed = 24;
    cellInfo(11).stimElec       = 23;
    cellInfo(11).type           = 'onMidg';

    %axonal stim well above threshold
    cellInfo(12).targetCell     = 348;
    cellInfo(12).otherCells     = 811;
    cellInfo(12).otherStimAxon  = true;
    %cellInfo(12).otherCells     = [766 811]; %766 doesn't show up on plot (above analyzed region)
    %cellInfo(12).otherStimAxon  = [true true];
    cellInfo(12).maxMovAnalyzed = 21;
    cellInfo(12).stimElec       = 24;
    cellInfo(12).type           = 'onMidg';

    %only gets to 0.28 in tested range (thresh = 1.6)
    cellInfo(13).targetCell     = 379;
    cellInfo(13).otherCells     = [];
    cellInfo(13).otherStimAxon  = [];
    cellInfo(13).maxMovAnalyzed = 26;
    cellInfo(13).stimElec       = 26;
    cellInfo(13).type           = 'onMidg';

    %axon stim below thresh
    cellInfo(14).targetCell     = 395;
    cellInfo(14).otherCells     = 317;
    cellInfo(14).otherStimAxon  = true;
    %cellInfo(14).otherCells     = [317 754]; %754 doesn't show up on plot (above analyzed region)
    %cellInfo(14).otherStimAxon  = [true true];
    cellInfo(14).maxMovAnalyzed = 25;
    cellInfo(14).stimElec       = 27;
    cellInfo(14).type           = 'onMidg';

    %no stimulation but on corner of array
    cellInfo(15).targetCell     = 437;
    cellInfo(15).otherCells     = [];
    cellInfo(15).otherStimAxon  = [];
    cellInfo(15).maxMovAnalyzed = [];
    cellInfo(15).stimElec       = [];
    cellInfo(15).type           = 'onMidg';

    %no evidence of stimulation but can't fully analyze on one electrode and on
    %edge of array
    cellInfo(16).targetCell     = 496;
    cellInfo(16).otherCells     = [];
    cellInfo(16).otherStimAxon  = [];
    cellInfo(16).maxMovAnalyzed = [];
    cellInfo(16).stimElec       = [];
    cellInfo(16).type           = 'onMidg';

    %gets to 0.66 before 2 other cells start to respond; analyzable up to m26
    %if cell 496 isn't considered
    cellInfo(17).targetCell     = 512;
    cellInfo(17).otherCells     = [513 676];
    cellInfo(17).otherStimAxon  = [false true];
    cellInfo(17).maxMovAnalyzed = 22;
    cellInfo(17).stimElec       = 35;
    cellInfo(17).type           = 'onMidg';

    %not stimulated, but can't fully analyze one of the stim electrodes
    cellInfo(18).targetCell     = 526;
    cellInfo(18).otherCells     = [];
    cellInfo(18).otherStimAxon  = [];
    cellInfo(18).maxMovAnalyzed = [];
    cellInfo(18).stimElec       = [];
    cellInfo(18).type           = 'onMidg';

    %not stimulated
    cellInfo(19).targetCell     = 527;
    cellInfo(19).otherCells     = [];
    cellInfo(19).otherStimAxon  = [];
    cellInfo(19).maxMovAnalyzed = [];
    cellInfo(19).stimElec       = [];
    cellInfo(19).type           = 'onMidg';

    %selective up to 0.56 (highest amplitude current)
    cellInfo(20).targetCell     = 541;
    cellInfo(20).otherCells     = [];
    cellInfo(20).otherStimAxon  = [];
    cellInfo(20).maxMovAnalyzed = 26;
    cellInfo(20).stimElec       = 37;
    cellInfo(20).type           = 'onMidg';

    %no evidence of stimulation but can't analyze past .45 (m17) on main
    %electrode and can't complete analysis on one of the neighboring electrodes
    %either
    cellInfo(21).targetCell     = 590;
    cellInfo(21).otherCells     = [];
    cellInfo(21).otherStimAxon  = [];
    cellInfo(21).maxMovAnalyzed = [];
    cellInfo(21).stimElec       = [];
    cellInfo(21).type           = 'onMidg';

    % only gets to ~0.35 by m24, can't analyze cell 590 beyond m24
    cellInfo(22).targetCell     = 616;
    cellInfo(22).otherCells     = [];
    cellInfo(22).otherStimAxon  = [];
    cellInfo(22).maxMovAnalyzed = 24;
    cellInfo(22).stimElec       = 42;
    cellInfo(22).type           = 'onMidg';

    % only gets to about 30%, can't fully analyze for cell 437
    cellInfo(23).targetCell     = 646;
    cellInfo(23).otherCells     = [];
    cellInfo(23).otherStimAxon  = [];
    cellInfo(23).maxMovAnalyzed = 26;
    cellInfo(23).stimElec       = 44;
    cellInfo(23).type           = 'onMidg';

    %full curve, fully selective
    cellInfo(24).targetCell     = 676;
    cellInfo(24).otherCells     = [];
    cellInfo(24).otherStimAxon  = [];
    cellInfo(24).maxMovAnalyzed = 14;
    cellInfo(24).stimElec       = 46;
    cellInfo(24).type           = 'onMidg';

    %unselective (figure cell); p49 slightly more selective but has higher
    %threshold
    cellInfo(25).targetCell     = 677;
    cellInfo(25).otherCells     = 317;
    cellInfo(25).otherStimAxon  = true;
    cellInfo(25).maxMovAnalyzed = 26;
    cellInfo(25).stimElec       = 46;
    cellInfo(25).type           = 'onMidg';

    %not stimulated, but cell on corner of array
    cellInfo(26).targetCell     = 707;
    cellInfo(26).otherCells     = [];
    cellInfo(26).otherStimAxon  = [];
    cellInfo(26).maxMovAnalyzed = 26;
    cellInfo(26).stimElec       = [];
    cellInfo(26).type           = 'onMidg';

    %mostly selective (small overlap)
    cellInfo(27).targetCell     = 766;
    cellInfo(27).otherCells     = 754;
    cellInfo(27).otherStimAxon  = false;
    cellInfo(27).maxMovAnalyzed = 21;
    cellInfo(27).stimElec       = 53;
    cellInfo(27).type           = 'onMidg';

    %selective up to 0.9
    cellInfo(28).targetCell     = 782;
    cellInfo(28).otherCells     = 858;
    cellInfo(28).otherStimAxon  = false;
    cellInfo(28).maxMovAnalyzed = 26;
    cellInfo(28).stimElec       = 54;
    cellInfo(28).type           = 'onMidg';

    %selective for full measured range (up to 90%)
    cellInfo(29).targetCell     = 811;
    cellInfo(29).otherCells     = [];
    cellInfo(29).otherStimAxon  = [];
    cellInfo(29).maxMovAnalyzed = 26;
    cellInfo(29).stimElec       = 55;
    cellInfo(29).type           = 'onMidg';

    %no evidence of stimulation, but can't confirm beyond m22 (0.73 µA) on one
    %elec and can't confirm in m26 for another elec
    cellInfo(30).targetCell     = 857;
    cellInfo(30).otherCells     = [];
    cellInfo(30).otherStimAxon  = [];
    cellInfo(30).maxMovAnalyzed = [];
    cellInfo(30).stimElec       = [];
    cellInfo(30).type           = 'onMidg';

    %partially selective - can't tell if other stim is axonal or somatic
    %(both cells very borderline)
    cellInfo(31).targetCell     = 858;
    cellInfo(31).otherCells     = [782 886];
    cellInfo(31).otherStimAxon  = [true false];
    cellInfo(31).maxMovAnalyzed = 26;
    cellInfo(31).stimElec       = 58;
    cellInfo(31).type           = 'onMidg';

    %fully selective
    cellInfo(32).targetCell     = 872;
    cellInfo(32).otherCells     = [];
    cellInfo(32).otherStimAxon  = [];
    cellInfo(32).maxMovAnalyzed = 26;
    cellInfo(32).stimElec       = 61;
    cellInfo(32).type           = 'onMidg';

    %mostly selective, axon stim at end
    cellInfo(33).targetCell     = 886;
    cellInfo(33).otherCells     = 348;
    cellInfo(33).otherStimAxon  = true;
    cellInfo(33).maxMovAnalyzed = 26;
    cellInfo(33).stimElec       = 60;
    cellInfo(33).type           = 'onMidg';

    %selective but on edge of mosaic
    cellInfo(34).targetCell     = 931;
    cellInfo(34).otherCells     = [];
    cellInfo(34).otherStimAxon  = [];
    cellInfo(34).maxMovAnalyzed = 26;
    cellInfo(34).stimElec       = 63;
    cellInfo(34).type           = 'onMidg';

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
for ii = 1:nCells
    load([pathToAnalysis 'elecResp_n' num2str(cellInfo(ii).targetCell) '_p' num2str(cellInfo(ii).stimElec)])
    
    elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', 0); %bootstrapReps = 100
    
    cellInfo(ii).thresh = elecResp.analysis.threshold;
    cellInfo(ii).sigma = elecResp.analysis.threshold/(-sqrt(2)*elecResp.analysis.erfParams(2));
    
    cellInfo(ii).maxAmpAnalyzed = abs(elecResp.stimInfo.stimAmps(elecResp.stimInfo.movieNos == cellInfo(ii).maxMovAnalyzed));
    
    %convert values to units of standard deviation relative to threshold
    cellInfo(ii).maxAnalyzedInSD = (cellInfo(ii).maxAmpAnalyzed - cellInfo(ii).thresh)/cellInfo(ii).sigma;
    
    %save([pathToAnalysis 'elecResp_n' num2str(cellInfo(ii).targetCell) '_p' num2str(cellInfo(ii).stimElec)], 'elecResp')
    
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
        load([pathToAnalysis 'elecResp_n' num2str(cellInfo(ii).otherCells(jj)) '_p' num2str(cellInfo(ii).stimElec)])
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
        
        %save([pathToAnalysis 'elecResp_n' num2str(cellInfo(ii).otherCells(jj)) '_p' num2str(cellInfo(ii).stimElec)], 'elecResp')
        
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

normInv30 = norminv(0.8, 0, 1); %distance from 0.5 to 0.8 probability, in SDs
normInv45 = norminv(0.95, 0, 1);%distance from 0.5 to 0.95 probability, in SDs

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


if ~plotVsCurrent
    % sort by maximum amplitude with verified selectivity
    maxSelective = [cellInfo.maxSelInSD];
    [~, plotOrder] = sort(maxSelective, 2, 'descend');
    
    %remove repeated fully selective cells %%%WARNING - if definition of fully
    %selective changes these numbers must change!! (but it should be obvious
    %from plot)
    plotOrder(1:4) = [];
    nPlottedRows = length(plotOrder);
else
    %order cells by target cell threshold
    targThreshs = [cellInfo.thresh];
    [~, plotOrder] = sort(targThreshs, 2, 'descend');
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
    %7 division because that's what's required to fit everything in 2011-01-11-0 selectivity summary
    barH = 0.8/7; %height of each non-target cell bar
    
    for ii = 1:nCells
        plotPos = find(plotOrder==ii);
        
        if ~isempty(plotPos)
            analyzedLim = cellInfo(ii).maxAnalyzedInSD;
            
            %plot bar
            plot([-normInv30 min([analyzedLim normInv30])], plotPos+1-0.5*barH*[1 1], 'k-', 'linewidth', 2)
            
            %plot circle for threshold
            if analyzedLim >= 0
                plot(0, plotPos+1-0.5*barH, 'o', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])
                %plot([0 0], plotPos+[1-0.02 1-barH+0.02], '-', 'lineWidth', 1, 'color', [0 0 0])
            end
            
            %plot triangle showing how far analysis was done
            if analyzedLim < normInv45
                plot(analyzedLim*[1 1], plotPos+[1 1], 'v', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])
            end
            
            %text(2, plotPos, ['n' num2str(cellInfo(ii).targetCell) ' p' num2str(cellInfo(ii).stimElec)])
            
            for jj = 1:length(cellInfo(ii).otherCells)
                
                curveRange = cellInfo(ii).otherThresh(jj) + cellInfo(ii).otherSigma(jj)*normInv30*[-1 1];
                
                %converts to units of target cell SDs from target cell thresh
                curveRange = (curveRange - cellInfo(ii).thresh)/cellInfo(ii).sigma;
                otherThresh = (cellInfo(ii).otherThresh(jj) - cellInfo(ii).thresh)/cellInfo(ii).sigma;
                
                threshColor = [1 0 0];
                
                yPlotVals = plotPos + [-0.02 -barH+0.02 -barH+0.02 -0.02] + 1 - jj*barH;
                
                %plot non-target cells
                if curveRange(2) <= -normInv45 %completely off left side of plot
                    plot(-normInv45-0.1, mean(yPlotVals), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', threshColor)
                elseif curveRange(1) <= -normInv45 %partially off left side of plot
                    if otherThresh < -normInv45
                        plot(-normInv45-0.1, mean(yPlotVals), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', threshColor)
                    end
                    plot([-normInv45-0.1 curveRange(2)], mean(yPlotVals)*[1 1], 'r-', 'linewidth', 2)
                elseif curveRange(2) <= analyzedLim %fully on plot
                    plot([curveRange(1) curveRange(2)], mean(yPlotVals)*[1 1], 'r-', 'linewidth', 2)
                    
                elseif curveRange(1) <= analyzedLim %partially off right side of plot
                    plot([curveRange(1) analyzedLim], mean(yPlotVals)*[1 1], 'r-', 'linewidth', 2)
                end
                
                if otherThresh <= analyzedLim && otherThresh >= -normInv45
                    plot(otherThresh, mean(yPlotVals), 'o', 'MarkerFaceColor', threshColor, 'MarkerEdgeColor', 'none')
                end
            end
        end
    end
else %plot in terms of current
    plotPosCurrent = 0;
    for ii = 1:length(plotOrder)
        cellInd = plotOrder(ii);
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
%         if otherCellLowerLim >= t+normInv30*sig && lim >= t+normInv30*sig %verified non-overlapping
%             plot(0.1, plotPosCurrent, 'k*')
%         elseif otherCellLowerLim >= t+normInv30*sig %unverified non-overlapping
%             plot(0.1, plotPosCurrent, 'k+')
%         end
        
        %plot(lim*[1 1], plotPosCurrent+[1 1], 'v', 'markerEdgeColor', 'none', 'markerFaceColor', [0 0 0])
        
        plotPosCurrent = plotPosCurrent+1;
    end
end

hold off

if ~plotVsCurrent
    %19 = nPlottedRows from 2011-01-11-0 selectivity summary--used here to
    %match row heights
    set(gca, 'xlim', [-normInv45-0.2 normInv45], 'ylim', [1 19+1.5], 'ytick', [], 'xtick', [-normInv45 -normInv30 0 normInv30 normInv45], 'xticklabel', [0.05 0.2 0.5 0.8 0.95])
    
    xlabel('target cell response probability')
else
    xlabel('current amplitude (µA)')
    set(gcf, 'color', [1 1 1], 'position', [100 0 400 8*plotPosCurrent+100])
    set(gca, 'ytick', [], 'ycolor', [1 1 1], 'ylim', [0 plotPosCurrent],...
        'xlim', [0 1.2], 'xtick', [0 0.4 0.8 1.2], 'units', 'pixels', 'position', [40 50 320 8*plotPosCurrent])
end

%% mosaic plots
datarun = load_data('2008-08-27-2/data001-lh/data001-lh');
datarun = load_index(datarun);
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));
datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
%datarun = get_sta_fits_from_vision(datarun,'all');
datarun.names.obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh';
datarun = load_obvius_sta_fits(datarun);
datarun = get_sta_fits_from_obvius(datarun, 'all');


% for electrode positions
datarun = compute_monitor_to_array_transformation(datarun);


%these lists should contain all cells were checked for nontarget activation, 
%which are all of the parasol and midget cells in vision analysis
%except those under "axon" or "off array" categories
goodCells{1} = [108 257 317 513 754 902]; %ON parasol
goodCells{2} = []; %OFF parasol
goodCells{3} = [18 91 92 138 166 167 182 241 256 274 332 348 379 395 496 512 ...
    526 527 541 590 616 646 676 677 707 766 782 811 857 858 872 886 931]; %ON midget
goodCells{4} = []; %OFF midget

stimElec = 8;
movieNo = 26;

selectivity_mosaic_plot(datarun, goodCells, stimElec, movieNo, pathToAnalysis)



%% plot of selectivity vs. threshold

PW = 100;
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
xlabel('threshold charge amplitude (pC)')
ylabel('selectivity index (target cell response prob.)')
set(gca,'xlim',[0 250], 'ylim', [0 1])


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

xlabel('threshold current amplitude (µA)')
ylabel('selectivity index (SDs from target cell thresh)')
set(gca,'xlim',[0 250], 'ylim', [-7 2])


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
































