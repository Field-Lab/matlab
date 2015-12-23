function cellTypeThreshCompareSummaryPlot(cellInfo, cellInfoNoStim, datarun, colors, examples, params, cellTypes)


% datarun must contain rf fits (used by plot_rf_summaries)
% colors must be 4x3 matrix specifying colors for each of the 4 cell types in this order: ON parasol, OFF parasol, ON midget, OFF midget
%
%
% params: structure with fields
%        .threshBarLim - lower and upper limits of threshold bar-plot (in
%        pC)  (default = [0 250])
%
%        .mosaicPlotLimX
%        .mosaicPlotLimY
%        .chargeLim - highest tested pulse amplitude, in pC
%        .safetyLim - location of 0.1 mC/cm^2, in pC
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('position', [100 100 800 700], 'color', [1 1 1])


%% mosaics plot

%calculate electrode positions
T_array_to_monitor = fliptform(datarun.piece.T_monitor_to_array);
T_monitor_to_STA = fliptform(coordinate_transform(datarun, 'monitor'));
T = maketform('composite',[T_monitor_to_STA T_array_to_monitor]);
e_positions = tformfwd(T, datarun.ei.position);

for ii = 1:length(cellTypes)
    cellsByType(ii).id = []; %#ok<AGROW>
end

for ii = 1:length(cellInfo)
    foundMatch = false;
    for jj = 1:length(cellTypes)
        if strcmpi(cellInfo(ii).type, cellTypes{jj})
            cellsByType(jj).id = [cellsByType(jj).id cellInfo(ii).id]; %#ok<AGROW>
            foundMatch = true;
            break
        end
    end
    if ~foundMatch
        error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
    end
end

for ii = 1:length(cellInfoNoStim)
    foundMatch = false;
    for jj = 1:length(cellTypes)
        if strcmpi(cellInfoNoStim(ii).type, cellTypes{jj})
            cellsByType(jj).id = [cellsByType(jj).id cellInfoNoStim(ii).id]; %#ok<AGROW>
            foundMatch = true;
            break
        end
    end
    if ~foundMatch
        error(['unrecognized cell type for neuron ' num2str(cellInfoNoStim(ii).id)])
    end
end

% plot on single array
if 0
    axes('position', [0.35 0.4 0.3 0.5], 'XLim', [8.3 22.2], 'YLim', [9.6, 23.4])
    hold on
    for ii = 1:length(cellTypes)
        plot_rf_summaries(datarun, [], 'array', true, 'array_color', [0 0 0])
        
        if isfield(params, 'mosaicLineWidths')
            plot_rf_summaries(datarun, [cellsByType(ii).id], 'fit_color', colors(ii,:), 'fit_width', params.mosaicLineWidths(ii))
        else
            plot_rf_summaries(datarun, [cellsByType(ii).id], 'fit_color', colors(ii,:), 'fit_width', 1.2)
        end
        plot(8.5, 24+2*ii, 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:), 'MarkerSize', 10)
        text(10, 24+2*ii, cellTypes{ii}, 'fontsize', 20)
    end
    
    hold off
    set(gca, 'XLim', params.mosaicPlotLimX, 'YLim', params.mosaicPlotLimY)
    axis equal
    axis off
end


% plot on multiple arrays
if 1
    axesPositions = {[0.30 0.68 0.20 0.24], [0.30 0.38 0.20 0.24], [0.5 0.68 0.20 0.24], [0.5 0.38 0.20 0.24]};
    for ii = 1:length(cellTypes)
        axes('position', axesPositions{ii}, 'XLim', [8.3 22.2], 'YLim', [9.6, 23.4])
        hold on
        plot(e_positions(:,1), e_positions(:,2), 'k.')

        plot_rf_summaries(datarun, [], 'array', true, 'array_color', [0 0 0])
        
        if isfield(params, 'mosaicLineWidths')
            plot_rf_summaries(datarun, [cellsByType(ii).id], 'fit_color', colors(ii,:), 'fit_width', params.mosaicLineWidths(ii))
        else
            plot_rf_summaries(datarun, [cellsByType(ii).id], 'fit_color', colors(ii,:), 'fit_width', 1.2)
        end
        %plot(8.5, 24+2*ii, 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:), 'MarkerSize', 10)
        %text(10, 24+2*ii, cellTypes{ii}, 'fontsize', 20)
                
        hold off
        set(gca, 'XLim', params.mosaicPlotLimX, 'YLim', params.mosaicPlotLimY)
        axis equal
        axis off
    end
end


%% response curves: examples from 4 cells

for ii = 1:4
    cellID = examples.ids(ii);
    
    iCell = find([cellInfo.id] == cellID);
        
    xProj = examples.xLims{ii}(1):0.01:examples.xLims{ii}(2);
    projection = 0.5 + 0.5*erf(cellInfo(iCell).params(1)*xProj + cellInfo(iCell).params(2));
    
    
    if ii == 1
        axes('position', [0.1 0.7 0.2 0.2], 'XLim', examples.xLims{1}, 'YLim', [0 1], 'box', 'on', 'Fontsize', 18)
    elseif ii == 2
        axes('position', [0.1 0.4 0.2 0.2], 'XLim', examples.xLims{2}, 'YLim', [0 1], 'box', 'on', 'Fontsize', 18)
    elseif ii == 3
        axes('position', [0.7 0.7 0.2 0.2], 'XLim', examples.xLims{3}, 'YLim', [0 1], 'box', 'on', 'Fontsize', 18)
    else
        axes('position', [0.7 0.4 0.2 0.2], 'XLim', examples.xLims{4}, 'YLim', [0 1], 'box', 'on', 'Fontsize', 18)
    end
    hold on
    
    foundMatch = false;
    for jj = 1:length(cellTypes)
        if strcmpi(cellInfo(iCell).type, cellTypes{jj})
            thisColor = colors(jj,:);
            foundMatch = true;
            break
        end
    end
    if ~foundMatch
        error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
    end
    plot(cellInfo(iCell).data(1,:), cellInfo(iCell).data(2,:), 'o', 'MarkerEdgeColor', thisColor, 'MarkerFaceColor', thisColor)
    plot(xProj, projection, 'LineWidth', 1, 'color', thisColor)
    plot(cellInfo(iCell).thresh*[1 1], [0 0.5], '--', 'color', thisColor)
    
    
    hold off
end


%% 1-D plot of thresholds, means and standard deviations

axes('position', [0.1 0.07 0.8 0.25])
make_1D_thresh_plot(cellInfo, colors, params, cellTypes)















