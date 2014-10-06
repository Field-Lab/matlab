%fits latencies to a curve and then plots a summary

%parameters for all datasets
clear all

%analysis parameters
params.recalcAll = false; %redo all response curve fits and standard deviation calculations
params.nBootstrapReps = 100;
params.maxSuccRateCutoff = 0.4;
params.removeOffArrayCells = true; params.cutoffMult = 0.5;
params.exclude30ArrayData = true;
params.includeDataset = [true true true true false true true];


params.cellTypes = {'onPar', 'offPar', 'onMidg', 'offMidg', 'sbc'};
%exampleCells = [2 18 676 470 768];
%exampleCells = [2 18 676 470 136];
exampleCells = [872 18 167 470 136];

targRespForLat = 0.5; %extract latencies from movie with response probability closest to this value 
binSize = 0.025; %in ms

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

%seed values for histogram fitting
t0Guess = 0.2;
tauGuess = 0.2;
n = 3;

excludeQuestionableCells = false;
pausePlots = true;
checkLatencyAccuracies = false; %displays plots of data traces overlaid on EI, similar to GUI, in order to verify that latencies are correct

%%% 2012-03-21: latencies were verified for all cells in datasets represented in
%%% cell_list_2008_08_27_2(), cell_list_2011_01_11_0(), cell_list_2011_10_25_4(), cell_list_2012_01_27_3() and
%%% cell_list_sbc(exclude30ArrayData=true), that satisfy removeSmallSigEdgeCells with cutoff = 0.5
%%% observation: in some movies (<20% of total), 1 or 2 spikes latencies
%%% are slightly off due to overlap with other spikes, and in some other
%%% movies (~15%?) signal is noisy and may cause errors in the +/-1 sample
%%% range


%% dataset details

cellInfo = load_analysis_all_paper1_cells(params);


%% get analysis

figure
axesH = axes;

unqualCount = 0;
for xx = 1:length(cellInfo)
    % load analysis from elecResp files
    for ii = 1:length(cellInfo{xx})

        load(cellInfo{xx}(ii).elecRespPath)
        
        %get movie numbers with of response rates closest to target value
        succRates = elecResp.analysis.successRates;
        
        succRateDiff = abs(succRates - targRespForLat);
        [sortedSuccDiff sortOrder] = sort(succRateDiff);
        iLatMov = sortOrder(1:2);
        %disp(num2str(sortedSuccDiff(1:3)))
        %disp('')
        
        cellInfo{xx}(ii).latMovSuccRates = succRates(iLatMov);
        
        %if excludeQuestionableCells == true, remove cells that require 2 templates for analysis or have shifted
        %upper template limit (both indications of potential errors in
        %latency estimates)
        exclude = false;
        if excludeQuestionableCells
            for jj = 1:length(iLatMov)
                template = squeeze(elecResp.cells.mainEI(elecResp.cells.recElec, :));
                templateMinPos = find(template == min(template));
                tempMinEnd = elecResp.analysis.details.tempOffsetWindow{iLatMov(jj)}(2) + templateMinPos;
                
                if ~isempty(elecResp.cells.active{iLatMov(jj)}) || tempMinEnd<35;
                    exclude = true;
                    break
                end
            end
        end
        
        
        if length(iLatMov) == 2 && ~exclude
            %visually inspect latency estimates to make sure they're OK --
            %this code was thrown together quickly so may not be robust to
            %all analysis situations, but should be obvious when it fails!
            if checkLatencyAccuracies
                for kk = 1:length(iLatMov)
                    if isempty(elecResp.cells.active{iLatMov(kk)}) %no other templates used in fitting
                        figure('position', [100 100 1200 500])
                        
                        latencies = elecResp.analysis.latencies{iLatMov(kk)};
                        dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo, elecResp.stimInfo.movieNos(iLatMov(kk)), 99999);
                        centerChannel = elecResp.cells.recElec;
                        
                        %plot of successes and failures
                        axes('position', [0.1 0.45 0.3 0.4])
                        title([elecResp.names.rrs_short_name ', movie ' num2str(elecResp.stimInfo.movieNos(iLatMov(kk)))])
                        subtractionVector = squeeze(mean(dataTraces(:, centerChannel, :), 1));
                        dataToPlot = zeros(size(dataTraces, 1), size(dataTraces, 3));
                        for i = 1:length(latencies)
                            dataToPlot(i, :) = squeeze(dataTraces(i, centerChannel, :)) - subtractionVector;
                        end
                        
                        hold on
                        for i = 1:length(latencies)
                            if ~any(latencies(i,:))
                                plot(dataToPlot(i, :), 'Color',[0.8 0.8 0.8])
                            end
                        end
                        for i = 1:length(latencies)
                            if any(latencies(i,:))
                                plot(dataToPlot(i, :), 'r')
                            end
                        end
                        
                        %stripped-down version of refreshEiPlots
                        failuresBin = (latencies == 0);
                        
                        unflaggedBin = elecResp.analysis.details.analysisFlags{iLatMov(kk)}==0;
                        failuresBinForArtifact = logical(failuresBin.*min(unflaggedBin,[],2));
                        
                        channelsToUse = getCluster(centerChannel);
                        nChannels = length(channelsToUse);
                        
                        eiData = elecResp.cells.mainEI;
                        eiYMin = 5000; eiYMax = -5000;
                        eiYMin = min([min(min(eiData(channelsToUse,:))), eiYMin]);
                        eiYMax = max([max(max(eiData(channelsToUse,:))), eiYMax]);
                        
                        ei = findStimEi(dataTraces, elecResp, elecResp.stimInfo.movieNos(iLatMov(kk)), 0, eiData, channelsToUse, true);
                        targetEi = eiData(channelsToUse, :);
                        targetEiMinPos = find(squeeze(targetEi(channelsToUse==centerChannel,:))...
                            ==min(squeeze(targetEi(channelsToUse==centerChannel,:)))); %position of minimum on primary electrode
                        
                        eiYMin = min([min(min(min(ei))), eiYMin]);
                        eiYMax = max([max(max(max(ei))), eiYMax]);
                        
                        for j = 1:size(ei,2)
                            axes('position', [(j-1)/(nChannels) 0.1 1/(nChannels+1) 0.3])
                            hold on
                            plot(squeeze(targetEi(j, targetEiMinPos-10:targetEiMinPos+15)), 'LineWidth', 2, 'color', [0.8 0 0])
                            for k = 1:length(latencies)
                                if latencies(k)
                                    plot(squeeze(ei(k,j,:)), 'color', [1 0 0])
                                end
                            end
                            hold off
                        end
                        pause
                        close(gcf)
                    else
                        warndlg(['need to check ' elecResp.names.rrs_short_name ', movie ' num2str(elecResp.stimInfo.movieNos(iLatMov(kk))) ' using full gui due to presence of other templates in fitting'])
                    end
                end
            end
            
            
            cellInfo{xx}(ii).latMovNo = elecResp.stimInfo.movieNos(iLatMov);
            
            %make sure corresponding movies are fully analyzed and finalized
            for jj = 1:length(iLatMov)
                if isempty(elecResp.analysis.type{iLatMov(jj)}) || ~elecResp.analysis.finalized(iLatMov(jj))
                    error(['movie ' num2str(elecResp.stimInfo.movieNos(iLatMov(jj))) ' from ' elecResp.names.rrs_short_name ' either analyzed or analysis isn''t finalized'])
                end
            end
            
            %subtract 2 from nonzero latency values because stimulus is applied at sample 2
            elecRespCopy = elecResp; %so that it isn't accidentally saved out with altered latency values
            for jj = 1:length(elecResp.stimInfo.movieNos)
                if ~isempty(elecResp.analysis.latencies{jj})
                    elecRespCopy.analysis.latencies{jj}(elecResp.analysis.latencies{jj}~=0) =...
                        elecResp.analysis.latencies{jj}(elecResp.analysis.latencies{jj}~=0)-2;
                end
            end
            
            cellInfo{xx}(ii).latencies = [];
            for jj = 1:length(iLatMov)
                cellInfo{xx}(ii).latencies = [cellInfo{xx}(ii).latencies; elecRespCopy.analysis.latencies{iLatMov(jj)}];
            end
            
            cellInfo{xx}(ii).nSpikes = sum(cellInfo{xx}(ii).latencies~=0);

                       
            histVals = plot_psth(axesH, elecRespCopy, cellInfo{xx}(ii).latMovNo, 'binSize', binSize);

            %shift time values so that they are in the center of
            %corresponding histogram bins            
            histVals(1,:) = histVals(1,:) + 0.5*(histVals(1,2) - histVals(1,1));
            
            %fit latencies to impulse response function
            [estParams alpha FWHM ttp fitError] = impRespFitter(histVals, t0Guess, tauGuess, n, 'makePlot', true, 'plotAxes', gca);
                       
            title([elecResp.names.rrs_short_name ', movie ' num2str(cellInfo{xx}(ii).latMovNo)...
                '(p=' num2str(sum(cellInfo{xx}(ii).latencies~=0)/length(cellInfo{xx}(ii).latencies)) ')'])
            set(gca, 'xlim', [0 2])
            
            if pausePlots
                pause
            end
            cla
            
            cellInfo{xx}(ii).estParams = [estParams alpha];
            cellInfo{xx}(ii).fitError = fitError;
            cellInfo{xx}(ii).ttp = ttp;
            cellInfo{xx}(ii).FWHM = FWHM;
        else
            unqualCount = unqualCount + 1;
            cellInfo{xx}(ii).latMovNo = [];
            cellInfo{xx}(ii).latencies = [];

            cellInfo{xx}(ii).estParams = [];
            if exclude
                disp(['*** excluded ' elecResp.names.rrs_short_name ' due to shifted template limit or presence of other neurons in analysis'])
            else
                disp(['*** excluded ' elecResp.names.rrs_short_name ' due to insufficient analysis'])
            end
        end
        
        clear elecResp succRates succDiff iMaxMov iLatMov
    end
end

close(gcf)


%% plot summary
cellTypes = {'on parasol', 'off parasol', 'on midget', 'off midget', 'sbc'};

tProj = 0:0.001:2;


figure('position', [100 100 400 400]); hold on;
errorSum = 0;
pooled_t0 = [];
pooled_tau = [];
pooled_latencies = [];
pooled_ttp = [];
pooled_FWHM = [];

for ii = 1:length(cellInfo)
    subplot(5,1,ii); hold on
    all_t0 = [];
    all_tau = [];
    all_latencies = [];
    all_ttp = [];
    all_FWHM = [];
    for jj = 1:length(cellInfo{ii})
        %plot histogram of example cells
        if cellInfo{ii}(jj).id == exampleCells(ii); %only works if there is only 1 cell with this ID within this type, but should be obvious when it doesn't work!
            binEdges = 0:binSize:2; %in ms
            psthPlotterBase(gca, binEdges, cellInfo{ii}(jj).latencies/20, 'normalizeToPeak', true, 'fillHist', true, 'lineColor', [0.7 0.7 0.7]);
            curveColor = [0 0 0];
            curveWidth = 2;
        else
            curveColor = [0.3 0.3 0.3];
            curveWidth = 0.5;
            %curveColor = [0.5 0.5 0.5];
        end
        
        %plot fit
        t0 = cellInfo{ii}(jj).estParams(1);
        tau = cellInfo{ii}(jj).estParams(2);
        t_sc = (tProj - t0)/tau; %scaled/shifted time values
        p = exp(-n*(t_sc-1)).*t_sc.^n;
        p(t_sc<0) = 0; %zero out values corresonding to t_xc < 0
        plot(tProj, p, 'color', curveColor, 'linewidth', curveWidth)
        %fill(tProj, p, colors(ii,:), 'edgeColor', 'none')
                     
        all_t0 = [all_t0 t0];
        all_tau = [all_tau tau];
        all_latencies = [all_latencies; cellInfo{ii}(jj).latencies];
        all_FWHM = [all_FWHM cellInfo{ii}(jj).FWHM];
        all_ttp = [all_ttp cellInfo{ii}(jj).ttp];
        errorSum = errorSum + cellInfo{ii}(jj).fitError;
    end
    %plot curve based on mean of each nonlinear parameter
%     t0 = mean(all_t0);
%     tau = mean(all_tau);
%     t_sc = (tProj - t0)/tau; %scaled/shifted time values
%     p = exp(-n*(t_sc-1)).*t_sc.^n;
%     p(t_sc<0) = 0; %zero out values corresonding to t_xc < 0
%     plot(tProj, p, 'color', [0 0 0], 'lineWidth', 2)
    ttp_mean = mean(all_ttp);
    ttp_SEM = std(all_ttp)/sqrt(length(all_ttp));
    
    lat_mean = mean(all_t0);
    lat_SEM = std(all_t0)/sqrt(length(all_t0));
    
    FWHM_mean = mean(all_FWHM);
    FWHM_SEM = std(all_FWHM)/sqrt(length(all_FWHM));

    %text(0,0.8, ['ttp = ' num2str(ttp_mean) '+/-' num2str(ttp_SEM)])
    %text(0,0.5, ['latency = ' num2str(lat_mean) '+/-' num2str(lat_SEM)])
    %text(0,0.2, ['FWHM = ' num2str(FWHM_mean) '+/-' num2str(FWHM_SEM)])

    title(cellTypes{ii})
    set(gca, 'ytick', [], 'ylim', [0 1.2], 'xlim', [0 1])
    
    pooled_t0 = [pooled_t0 all_t0];
    pooled_ttp = [pooled_ttp all_ttp];
    pooled_FWHM = [pooled_FWHM all_FWHM];
end
ylim = get(gca, 'ylim');
%text(1,0.8*ylim(2), ['error = ' num2str(errorSum)])


pooled_t0_mean =   mean(pooled_t0);
pooled_t0_SEM =    std(pooled_t0)/sqrt(length(pooled_t0));
pooled_ttp_mean =  mean(pooled_ttp);
pooled_ttp_SEM =   std(pooled_ttp)/sqrt(length(pooled_ttp));
pooled_FWHM_mean = mean(pooled_FWHM);
pooled_FWHM_SEM =  std(pooled_FWHM)/sqrt(length(pooled_FWHM));





