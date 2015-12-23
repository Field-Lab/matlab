%% collection of plotting scripts for primate cluster stimulation summaries

clear all

reExtractAnalysis = false;

%plotting parameters
plotParams.showSecAloneThresh = true;
plotParams.checkForPAloneArt = false;

useFullModel = true; %whether to calculate estimated shifts from full model (using lambdas), vs. measured component pair shifts
useFitPAloneThresh = true; %whether to use primary-alone threshold estimated from linear fits of pair pattern thresholds rather than actual measured primary-alone threshold
plotUnfittedTriplets = false; %whether to plot the triplets corresponding to unfitted component pairs in obs vs expected plot

%fitting parameters
constrainSlopes = 0;
plotInconsistentFits = false;
params.useRobustFitting = true;
params.refitAllErf = false;

linErrThresh = 0.0002; %threshold error value for pair data for a particular secondary electrode to be split into multiple regions

params.requiredSuccessRate =  .4; % a "minimum possible threshold" is calculated if there are no measured success rates above this value (plotted in addition to secondary-alone threhsolds)

params.checkForAnalysisGaps = false;

badElecThresh = 0.5; %fraction of mean(over electrodes) summed EI amplitude that an electrode
% must fall below to be considered "bad"

suppressIndPlots = false; %suppress plotting of linearity results for each individual cell

%%

%load list of cells
CI = cell_list_cluster_stim_StanfordDirs(badElecThresh);
nCells = length(CI);

%% loading/initial fitting of pairs data

for ii = 1:nCells
    disp(['cell ' num2str(ii)])
    % get ei from elecResp file
    tmp = load([CI(ii).pathToData filesep 'elecResp_n' num2str(CI(ii).id) '_p' num2str(CI(ii).patternNos(1))]);
    CI(ii).ei = tmp.elecResp.cells.mainEI;    
    
    % extracts data 
    clusterFN = [CI(ii).pathToData filesep 'n' num2str(CI(ii).id) '_clusterData.mat'];
    if ~exist(clusterFN, 'file') || reExtractAnalysis
        [CI(ii).thresholds, CI(ii).threshStds, CI(ii).details] = extractPatternStimAnalysis(CI(ii).pathToData, CI(ii).patternNos,...
            CI(ii).pElec, CI(ii).id, 'excludeFromFits', CI(ii).excludePatterns, params);
        
        % store extracted data for later use
        clusterData.thresholds = CI(ii).thresholds;
        clusterData.threshStds = CI(ii).threshStds;
        clusterData.details = CI(ii).details;
        
        save(clusterFN, 'clusterData')
    else
        load(clusterFN)
        CI(ii).thresholds = clusterData.thresholds;
        CI(ii).threshStds = clusterData.threshStds;
        CI(ii).details = clusterData.details;
    end
        
    % fit pair data
    [CI(ii).details.tp CI(ii).details.lambdas CI(ii).details.unfitPairs CI(ii).details.mValidBound CI(ii).details.multiRegModel] =...
        linFitSimulMinError(CI(ii).thresholds, CI(ii).details, 'linErrThresh', linErrThresh);
        
    % save first-order analysis
    CIFN = [CI(ii).pathToData filesep 'n' num2str(CI(ii).id) '_CIData.mat'];
    CIData = CI(ii);
    %save(CIFN, 'CIData')
end

%% tool to test value of linErrThresh and make histogram figure of nonlinearity indeces
%****** ONLY WORKS FOR ORIGINAL 15 CELLS (hand-labeled)

checkLinErrMeasures(CI, linErrThresh)


%% load up saved CI array with fits
%load('Users/lhruby/Desktop/CItmp.mat')

            
%% plotting/triplet analysis

for ii = 1:nCells
    if ~suppressIndPlots
        hAxes = makeClusterStimEILayout(CI(ii).ei, CI(ii).pElec, 'neuronID', CI(ii).id, 'pathToData', CI(ii).pathToData, 'cellType', CI(ii).type);

        pairPatternPlots(CI(ii).thresholds, CI(ii).threshStds, CI(ii).details, 'xPlotLim', CI(ii).xPlotLim, 'yPlotLim', CI(ii).yPlotLim,...
            'plotAxes', hAxes, 'plotUnfitPairs', false, plotParams);
    end
                
    % need to make "suppressIndPlots" a param-value for
    % generateSecondaryCombPlots
    if CI(ii).tripletsAnalyzed && ~suppressIndPlots
        figure('position', [100 100 500 500])
        hAxes = axes('position', [0.05 0.05 0.9 0.9]); hold on
        keyboard;
         [r_sq3 CI(ii).details.unfitTrip] = generateSecondaryCombPlots(CI(ii).thresholds, CI(ii).threshStds, CI(ii).details, 'fullModel', true, 'useFitPAloneThresh', true, 'obsVExpAxes', hAxes, 'plotUnfitted', plotUnfittedTriplets);
%         generateSecondaryCombPlots(secondaryCombs, thresholds, threshStds, estPrimThresh, actualRelAmps)
         %r_sq1 = generateSecondaryCombPlots(CI(ii).thresholds, CI(ii).threshStds, CI(ii).details, 'fullModel', false, 'useFitPAloneThresh', false, 'obsVExpAxes', hAxes, 'plotUnfitted', plotUnfittedTriplets, 'obsVExpColor', [1 0 0]);
         %r_sq2 = generateSecondaryCombPlots(CI(ii).thresholds, CI(ii).threshStds, CI(ii).details, 'fullModel', true, 'useFitPAloneThresh', false, 'obsVExpAxes', hAxes, 'plotUnfitted', plotUnfittedTriplets, 'obsVExpColor', [0 0 1]);
         %[r_sq3 CI(ii).details.unfitTrip] = generateSecondaryCombPlots(CI(ii).thresholds, CI(ii).threshStds, CI(ii).details, 'fullModel', true, 'useFitPAloneThresh', true, 'obsVExpAxes', hAxes, 'plotUnfitted', plotUnfittedTriplets, 'obsVExpColor', [0 1 0]);
% 
         %title(['full model, measured p-alone thresh (blue, R^2 = ' num2str(r_sq2)  ')' 10 'vs full model, p-alone thresh from fits (black, R^2 = ' num2str(r_sq3) ') vs additivity only (red, R^2 = ' num2str(r_sq1)  ')'])
%         xl = get(hAxes, 'xlim'); yl = get(hAxes, 'ylim');
%         text(xl(1)+(xl(2)-xl(1))*0.05, yl(2)-(yl(2)-yl(1))*0.05, ['neuron', num2str(CI(ii).id), 10, CI(ii).pathToData])
%         text(xl(1)+(xl(2)-xl(1))*0.05, yl(2)-(yl(2)-yl(1))*0.15,...
%             ['full model, measured p-alone thresh (blue, R^2 = ' num2str(r_sq2)  ')' 10,...
%             'vs full model, p-alone thresh from fits (green, R^2 = ' num2str(r_sq3) ')', 10,...
%             'vs additivity only (red, R^2 = ' num2str(r_sq1)  ')'])

    elseif CI(ii).tripletsAnalyzed
        figure('position', [100 100 250 250])
        hAxes = axes('position', [0.2 0.2 0.6 0.6]); hold on
        [r_sq, CI(ii).details.unfitTrip] = generateSecondaryCombPlots(CI(ii).thresholds, CI(ii).threshStds, CI(ii).details,...
            'fullModel', useFullModel, 'useFitPAloneThresh', useFitPAloneThresh, 'obsVExpAxes', hAxes, 'plotUnfitted', plotUnfittedTriplets, 'obsVExpColor', [1 0 0]);
        title(['cell ' num2str(CI(ii).id), 10, 'R^2 = ' num2str(r_sq)])
    end
end


%% second-order analysis plots

%just to get .meanAxAng field
CI = plotLambdaVAxonDir(CI, 'excludeBadElecs', true, 'range90', false, 'binSize', 30);

CI_60 = CI([CI.pitch]==60);
CI_30 = CI([CI.pitch]==30);

plotLambdaVAxonDir(CI_60, 'excludeBadElecs', false, 'binSize', 60, 'excludeCellsInd', 4, 'nShuffleReps', 50000, 'lamLimPlot', 0.5);
plotLambdaVAxonDir(CI_30, 'excludeBadElecs', false, 'binSize', 60, 'nShuffleReps', 50000);

%% slope vs array spacing
plotSlopeVArraySpace(CI, 'excludeBadElecs', true, 'binSpacing', 0.1)

%% slope vs EI signal

%plotLambdaVEISig(CI_60, 'useAbsEISig', false, 'nShuffleReps', 1000, 'excludeCellsInd', 4, 'lambdaType', 'positive')
%plotLambdaVEISig(CI_60, 'useAbsEISig', false, 'nShuffleReps', 1000, 'excludeCellsInd', 4, 'lambdaType', 'negative')
plotLambdaVEISig(CI_60, 'useAbsEISig', false, 'nShuffleReps', 1000, 'excludeCellsInd', 4, 'lambdaType', 'all', 'excludeBadElecs', false)

%only use 'all' because there are no negative lambda values
plotLambdaVEISig(CI_30, 'useAbsEISig', false, 'nShuffleReps', 1000, 'lambdaType', 'all', 'excludeBadElecs', false)

plotLambdaVEISig(CI, 'useAbsEISig', false, 'nShuffleReps', 1000, 'excludeCellsInd', 5, 'lambdaType', 'all', 'excludeBadElecs', false)


%% special plots for selectivity improvement example
% FIGURE 6
selImpExp = '2011-01-11-0/data032';
paths = {CI.pathToData};

dsMatch = cellfun(@(x) ~isempty(strfind(x, selImpExp)), paths);

bestPInd = [102 118 190 120 192];
plotSelectivityImprovementExample(CI(dsMatch), bestPInd)

%% plot lambdas over ei movie to get a sense for how lambda values relate
% to axon orientation (saves movie)

for ii = 1:nCells
%ii = 1
    mov = plotLambdaOnEIMovie(CI(ii));
    save(['/Users/lhruby/Desktop/eiMovies/eiMovie' num2str(CI(ii).id)], 'mov')
    clear mov
end
%% play the movies you made

figure('position', [100 100 800 800])

for ii = 1
    load(['/Users/lhruby/Desktop/eiMovies/eiMovie' num2str(CI(ii).id)])
    for jj = 1:1
        movie(mov, 1, 15)
    end
    clear mov
end



%% example EIs/lambdas showing how lambda values are related to axon orientation

CIToPlot = [1 6 11 13 3];
[xCoords yCoords] = getElectrodeCoords61();
scaleFactorLam = 75;

for ii = 1:length(CIToPlot)
    thisEI = CI(CIToPlot(ii)).ei;

    % calculate maximum waveform value on each electrode (absolute value)
    eiAmps = zeros(64, 1);
    for j = 1:64
        if ~(j==9||j==25||j==57)
            eiAmps(j) = max(max(abs(thisEI(j,:))));
        end
    end
    eiAmps = eiAmps/max(eiAmps);
    
    contours = hex_contour(xCoords, yCoords, eiAmps, 12, 'fig_or_axes', 0, 'contourSpacing', 'linear');
    
    %array outline
    hold on
    plot([0 8.6603 8.6603 0 -8.6603 -8.6603 0], [10 5 -5 -10 -5 5 10], 'k-')
    
    %axon orientation and lambdas
    
    xAx = xCoords(CI(CIToPlot(ii)).pElec) + [0 cosd(CI(CIToPlot(ii)).meanAxAng)];
    yAx = yCoords(CI(CIToPlot(ii)).pElec) + [0 sind(CI(CIToPlot(ii)).meanAxAng)];
        
    for jj = 1:6
        sElec = CI(CIToPlot(ii)).details.sElecs(jj);
        if CI(CIToPlot(ii)).details.lambdas(jj) > 0
            plot(xCoords(sElec), yCoords(sElec), 'ko', 'markerSize', ceil(scaleFactorLam*CI(CIToPlot(ii)).details.lambdas(jj)), 'linewidth', 2)
        else
            plot(xCoords(sElec), yCoords(sElec), 'o', 'markerSize', ceil(scaleFactorLam*abs(CI(CIToPlot(ii)).details.lambdas(jj))),...
                'markerEdgeColor', [0.5 0.5 0.5], 'linewidth', 2)
        end
    end
    pElec = CI(CIToPlot(ii)).pElec;
    plot(xCoords(pElec), yCoords(pElec), 'k*', 'markerSize', 10)
    plot(xAx, yAx, 'k-', 'linewidth', 2)

    
    set(gca, 'XLim', [-10 10], 'YLim', [-11 11])
    title(['neuron' num2str(CI(CIToPlot(ii)).id)])
    axis equal
    axis off
end


