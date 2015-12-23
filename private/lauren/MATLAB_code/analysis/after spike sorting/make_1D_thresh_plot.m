function make_1D_thresh_plot(cellInfo, colors, params, cellTypes, markedCells)

if ~exist('markedCells', 'var')
    markedCells = [];
end

yBarPos = params.yBarPos;
OS = params.OS;

xLength = params.threshBarLim(2) - params.threshBarLim(1);

%gather up thresholds (in charge amplitude) for different groups
for ii = 1:length(cellTypes)
    thresh(ii).ver = []; %#ok<*AGROW>
    thresh(ii).unver = [];
    thresh(ii).marked = [];
end

for ii = 1:length(cellInfo)
    %determine pulse phase width
    if isfield(cellInfo(ii), 'PW') %use PW from cellInfo struct if specified
        PW = cellInfo(ii).PW;
        %make sure this is consistent with params.PW
        if length(params.PW) == 1 && params.PW ~= PW
            error('param.PW is not consistent with PW specified in cellInfo structure')
        end
    elseif length(params.PW) == 1 %otherwise use PW specified in params structure
        PW = params.PW;
    else %can't unambiguously determine PW
        error('param.PW is not a single value and PW is not specified in cellInfo structure')
    end
       
    foundMatch = false;
    for jj = 1:length(cellTypes)
        if strcmpi(cellInfo(ii).type, cellTypes{jj})
            if any(markedCells == cellInfo(ii).id)
                thresh(jj).marked = [thresh(jj).marked cellInfo(ii).thresh*PW];
            elseif cellInfo(ii).verMin
                thresh(jj).ver = [thresh(jj).ver cellInfo(ii).thresh*PW];
            else
                thresh(jj).unver = [thresh(jj).unver cellInfo(ii).thresh*PW];
            end
            foundMatch = true;
            break
        end
    end
    if ~foundMatch
        error(['unrecognized cell type for neuron ' num2str(cellInfo(ii).id)])
    end
end

for ii = 1:length(cellTypes)
    thresh(ii).all = [thresh(ii).ver thresh(ii).unver];
end


%axes('position', [0.1 0.07 0.8 0.25])

%plot bar
hold on
plot([params.threshBarLim(1) params.threshBarLim(2)], [0 0], 'k-', 'Linewidth', 1)
plot([params.threshBarLim(1) params.threshBarLim(1)], [-0.1 0.1], 'k-', 'Linewidth', 1)
plot([params.threshBarLim(2) params.threshBarLim(2)], [-0.1 0.1], 'k-', 'Linewidth', 1)
plot(0.5*[params.threshBarLim(2) params.threshBarLim(2)], [-0.1 0.1], 'k-', 'Linewidth', 1)





%plot individual thresholds
for ii=1:length(cellTypes)
    
    %plot grey line
    plot([params.threshBarLim(1) params.threshBarLim(2)], yBarPos(ii)*[1 1], 'color', [0.8 0.8 0.8], 'Linewidth', 1)
    
    %subcategory means, square style
    %plot(mean(thresh(ii).ver),   yBarPos(ii), 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10)
    %plot(mean(thresh(ii).unver), yBarPos(ii), 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 1], 'MarkerSize', 10)
    
    %plot individual thresholds
    plot(thresh(ii).ver,   zeros(size(thresh(ii).ver))+yBarPos(ii),   'o', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', 'none')
    %plot(thresh(ii).ver,   zeros(size(thresh(ii).ver))+yBarPos(ii),   'o', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:))
    plot(thresh(ii).unver, zeros(size(thresh(ii).unver))+yBarPos(ii), 'o', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', 'none')
    %plot(thresh(ii).unver, zeros(size(thresh(ii).unver))+yBarPos(ii), 'o', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:))
    
    plot(thresh(ii).marked, zeros(size(thresh(ii).marked))+yBarPos(ii), 'x', 'MarkerEdgeColor', colors(ii,:))
    
    
    %subcategory means, vertical bar style
    %plot(mean(thresh(ii).ver)*[1 1],   yBarPos(ii)+[-OS OS], '-', 'color', [0 0 0], 'Linewidth', 2)
    %plot(mean(thresh(ii).unver)*[1 1 1 1 1] + 0.002*xLength*[-1 -1 1 1 -1], yBarPos(ii)+[-OS OS OS -OS -OS], '-', 'color', [0 0 0], 'Linewidth', 1)
    
    
%     plot(mean(thresh(ii).ver),   yBarPos(ii)+OS, 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 8)
%     plot(mean(thresh(ii).unver), yBarPos(ii)+OS, 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'none', 'MarkerSize', 8)

%     %subcategory SDs
%     plot([mean(thresh(ii).ver) - std(thresh(ii).ver),  mean(thresh(ii).ver) + std(thresh(ii).ver)],...
%        (yBarPos(ii)+2*OS)*[1 1], 'color', colors(ii,:), 'linewidth', 1);
%     
%     plot([mean(thresh(ii).unver) - std(thresh(ii).unver),  mean(thresh(ii).unver) + std(thresh(ii).unver)],...
%        (yBarPos(ii)+OS)*[1 1], 'color', colors(ii,:), 'linewidth', 1);
%     
%     
%     
    
%     plot(thresh(ii).ver,   zeros(size(thresh(ii).ver)),   'o', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:))
%     plot(thresh(ii).unver, zeros(size(thresh(ii).unver)), 'o', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', 'none')
    
%     %plot means
%     plot(mean(thresh(ii).all),  yBarPos(ii), 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:), 'MarkerSize', 10)
% 
%     %standard deviations
%     plot([mean(thresh(ii).all) - std(thresh(ii).all),  mean(thresh(ii).all) + std(thresh(ii).all)],...
%         [yBarPos(ii) yBarPos(ii)], 'color', colors(ii,:), 'linewidth', 2);

    
%     %subcategories
%     plot(mean(thresh(ii).ver),  yBarPos(ii)+OS, 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', colors(ii,:), 'MarkerSize', 5)
%     plot([mean(thresh(ii).ver) - std(thresh(ii).ver),  mean(thresh(ii).ver) + std(thresh(ii).ver)],...
%         [yBarPos(ii)+OS yBarPos(ii)+OS], 'color', colors(ii,:), 'linewidth', 1);
%     
%     plot(mean(thresh(ii).unver),  yBarPos(ii)+2*OS, 's', 'MarkerEdgeColor', colors(ii,:), 'MarkerFaceColor', 'none', 'MarkerSize', 5)
%     plot([mean(thresh(ii).unver) - std(thresh(ii).unver),  mean(thresh(ii).unver) + std(thresh(ii).unver)],...
%         [yBarPos(ii)+2*OS yBarPos(ii)+2*OS], 'color', colors(ii,:), 'linewidth', 1);
    
end

%plot markers for tested current range and safety limitations
if isfield(params, 'chargeLim')
    if params.chargeLim <= params.threshBarLim(2)
        plot(params.chargeLim, max(yBarPos)+0.25, 'v', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none')
    end
end

if isfield(params, 'safetyLim')
    plot(params.safetyLim, 0, '^', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none')
end

if length(params.PW) == 1
    text(0, -0.5, '0 (0)', 'HorizontalAlignment', 'center', 'FontSize', 18)
    text(0.5*params.threshBarLim(2), -0.5, [num2str(0.5*params.threshBarLim(2)/params.PW) '(' num2str(0.5*params.threshBarLim(2)) ')'], 'HorizontalAlignment', 'center', 'FontSize', 18);
    text(params.threshBarLim(2), -0.5, [num2str(params.threshBarLim(2)/params.PW) '(' num2str(params.threshBarLim(2)) ')'], 'HorizontalAlignment', 'center', 'FontSize', 18);
else
    text(0, -0.5, '(0)', 'HorizontalAlignment', 'center', 'FontSize', 18)
    text(0.5*params.threshBarLim(2), -0.5, ['(' num2str(0.5*params.threshBarLim(2)) ')'], 'HorizontalAlignment', 'center', 'FontSize', 18);
    text(params.threshBarLim(2), -0.5, ['(' num2str(params.threshBarLim(2)) ')'], 'HorizontalAlignment', 'center', 'FontSize', 18);
end
text(0.5*params.threshBarLim(2), -1, 'threshold: µA (pC)', 'HorizontalAlignment', 'center', 'FontSize', 20);
hold off
set(gca, 'YLim', [-1 1.5], 'XLim', [params.threshBarLim(1) params.threshBarLim(2)])
axis off



