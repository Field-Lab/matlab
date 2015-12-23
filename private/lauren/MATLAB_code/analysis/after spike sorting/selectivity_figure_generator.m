%function selectivity_figure_generator(pathToElecResp, pathToVisResp, cellIDs, targetCell, patternNo, movieNo, plotElec)


%% SEE PLOTSELECTIVITY.M FOR NEWER VERSION



%% for testing
clear all

pathToElecResp = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002';
pathToVisResp = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh';

%pathToElecResp = '/Ripple/Analysis/lauren/2008-08-27-2/data002';
%pathToVisResp = '/Ripple/Analysis/lauren/2008-08-27-2/data001-lh/data001-lh';

cellIDs{1} = [18 91 92 138 166 167 182 241 256 274 332 348 379 395 ...
    496 512 526 527 541 590 616 646 676 677 707 766 782 811 857 858 872 886 931]; %437 left out due to insufficient signal size

cellIDs{2} = [108 257 317 513 754 902];

%targetCell = 677;
%targetCell = 167;

%patternNo = 46;
%movieNo = 23;
%patternNo = 49;
%movieNo = 26;
%patternNo = 12;
%movieNo = 19;

patternNo = 19;
movieNo = 19;
targetCell = 274;


plotElec = [];

%% load up visual responses

%datarun = load_data('2008-08-27-2/data001-lh/data001-lh');
datarun = load_data(pathToVisResp);


datarun = load_index(datarun);
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));
datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
datarun = load_ei(datarun,[]);
datarun = get_sta_fits_from_vision(datarun,'all');

% compute the transformation using the clicked points in datarun.piece.array
datarun = compute_monitor_to_array_transformation(datarun);

%% load up electrical responses

respProb = cell(size(cellIDs));

for jj = 1:length(cellIDs)
    respProb{jj} = zeros(size(cellIDs{jj}));
    for ii = 1:length(cellIDs{jj})
        temp = load([pathToElecResp '/elecResp_n' num2str(cellIDs{jj}(ii)) '_p' num2str(patternNo) '.mat']);
        elecResp = temp.elecResp;

        movieI = find(elecResp.stimInfo.movieNos == movieNo);
        
        %makes sure nPulses is correct
        if elecResp.stimInfo.nPulses(movieI) == 0 %value hasn't been set yet
            dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo, movieNo);
            elecResp.stimInfo.nPulses(movieI) = size(dataTraces, 1);
            elecResp.analysis.erfCurrent = 0;
        end
        
        %makes sure movie is analyzed and finalized, and successRate is correct
        if ~isempty(elecResp.analysis.type{movieI}) %movie has been analyzed
            successRate = sum(elecResp.analysis.latencies{movieI}~=0)/elecResp.stimInfo.nPulses(movieI);
            if elecResp.analysis.successRates(movieI) ~= successRate
                disp(['warning: success rate for movie ' num2str(movieNo) ' had to be corrected'])
                elecResp.analysis.successRates(movieI) = successRate;
                elecResp.analysis.erfCurrent = 0;
            end
            if ~elecResp.analysis.finalized(movieI)
                disp(['pattern ' num2str(patternNo) ', movie ' num2str(movieNo) ' has not been finalized for cell ' num2str(cellIDs{jj}(ii))])
            end
        else
            error(['pattern ' num2str(patternNo) ', movie ' num2str(movieNo) ' has not been analyzed for cell ' num2str(cellIDs{jj}(ii))])
        end

        respProb{jj}(ii) = elecResp.analysis.successRates(movieI);
    end
end


%% generate figure
colors = jet(21);
colors = colors(1:end-1, :); %20 colors that don't go into dark red

figure
for jj = 1:length(cellIDs)
    hAxes{jj} = axes('position', [0.05 + (jj-1)/length(cellIDs), 0.05, 1/(length(cellIDs)+0.5), 0.9]);
    hold on
    for ii = 1:20
        if ii == 1 %special case for first interval, since it includes the lower bound
            cellInds = find(respProb{jj} <= 1/20);
        else
            cellInds = find((respProb{jj} > (ii-1)/20) & (respProb{jj} <= ii/20));
        end
        plot_rf_summaries(datarun, cellIDs{jj}(cellInds), 'fit_width', 1, 'fit_color', colors(ii,:), 'array', false)
    end
    if ~any(cellIDs{jj}==targetCell) %target cell isn't in the current set (mosaic)
        plot_rf_summaries(datarun, targetCell, 'fit_width', 1, 'fit_color', [0.5 0.5 0.5], 'array', false)
    end
    
    hold off
    axis off
end


%plot electrode location if desired
for ii = 1:length(plotElec)
    for jj = 1:length(cellIDs)
        axes(hAxes{jj}); hold on
        plot(datarun.ei.position_sta(plotElec(ii), 1), datarun.ei.position_sta(plotElec(ii), 2), 'k.')
        hold off
    end
end


%plot vertical color scale bar

figure('position', [100 100 100 500])
hold on
nColors = size(colors, 1);
xPoly = [0 1 1 0];
for ii = 1:nColors
    yPoly = [(ii-1)/nColors (ii-1)/nColors (ii+0.05)/nColors (ii+0.05)/nColors];
    fill(xPoly, yPoly, colors(ii,:), 'EdgeColor', colors(ii,:), 'LineWidth', 0.001)
end

colorBarLabels = 0:0.5:1;

for ii = 1:length(colorBarLabels)
    y = (colorBarLabels(ii)-0)/(1-0);
    plot([1 1.1], [y y], 'k-')
    text(1.2, y, num2str(colorBarLabels(ii)))
end
set(gca, 'ylim', [0 1], 'xlim', [0 2])
axis off
hold off



