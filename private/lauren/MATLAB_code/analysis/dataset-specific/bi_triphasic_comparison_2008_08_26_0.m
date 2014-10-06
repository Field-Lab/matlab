% script for comparing biphasic and triphasic pulse responses in 2008-08-26-0/data002, data003
%function bi_triphasic_comparison_2008_08_26_0() %just so that there can be a subfunction

%%

clear all
cd /snle/lab/Experiments/Array/Analysis/2008-08-26-0/

%onParasolFiles = {'n948_p2', 'n137_p10', 'n272_p16', 'n303_p16', 'n303_p18', 'n272_p19', 'n558_p36',...
%    'n588_p39', 'n588_p40', 'n558_p44', 'n558_p47', 'n798_p54', 'n813_p55', 'n798_p56'};
onParasolFiles = {'n948_p2', 'n272_p16', 'n303_p16', 'n303_p18', 'n272_p19',...
    'n588_p40', 'n558_p44', 'n813_p55'};

%onParasolFiles1 = {'n558_p36', 'n558_p44'};
%onParasolFiles2 = {'n272_p16', 'n272_p19'};
%onParasolFiles3 = {'n303_p16', 'n303_p18'};
%onParasolFiles4 = {'n588_p40'};

offParasolFiles = {'n138_p10','n197_p14', 'n335_p23', 'n678_p46', 'n678_p49'};
%offParasolFiles1 = {'n678_p46', 'n678_p49'};

onMidgetFiles = {'n31_p3', 'n151_p7'};

%axonFiles = {'n497_p6','n497_p26', 'n678_p53'};
%axonFiles1 = {'n497_p6','n497_p26'};

sbcFiles = {'n466_p28'};

%% package filenames
allFiles = {};
fileCounter = 0;
cellIDCounter = 0;
cellType = [];
sameCells = [];

for ii = 1:length(onParasolFiles)
    fileCounter = fileCounter + 1;
    allFiles{fileCounter} = onParasolFiles{ii};
    cellType = [cellType 1];
    matchFound = false;
    for jj = 1:ii-1
        if strcmp(onParasolFiles{ii}(1:4), onParasolFiles{jj}(1:4))
            matchFound = true;
            sameCells = [sameCells sameCells(jj)];
        end
    end
    if ~matchFound
        cellIDCounter = cellIDCounter + 1;
        sameCells = [sameCells cellIDCounter];
    end
end

fileCountBtwnTypes = fileCounter;
for ii = 1:length(offParasolFiles)
    fileCounter = fileCounter + 1;
    allFiles{fileCounter} = offParasolFiles{ii};
    cellType = [cellType 2];
    matchFound = false;
    for jj = 1:ii-1
        if strcmp(offParasolFiles{ii}(1:4), offParasolFiles{jj}(1:4))
            matchFound = true;
            sameCells = [sameCells sameCells(jj+fileCountBtwnTypes)];
        end
    end
    if ~matchFound
        cellIDCounter = cellIDCounter + 1;
        sameCells = [sameCells cellIDCounter];
    end
end

fileCountBtwnTypes = fileCounter;
for ii = 1:length(onMidgetFiles)
    fileCounter = fileCounter + 1;
    allFiles{fileCounter} = onMidgetFiles{ii};
    cellType = [cellType 3];
    matchFound = false;
    for jj = 1:ii-1
        if strcmp(onMidgetFiles{ii}(1:4), onMidgetFiles{jj}(1:4))
            matchFound = true;
            sameCells = [sameCells sameCells(jj+fileCountBtwnTypes)];
        end
    end
    if ~matchFound
        cellIDCounter = cellIDCounter + 1;
        sameCells = [sameCells cellIDCounter];
    end
end

fileCountBtwnTypes = fileCounter;
for ii = 1:length(sbcFiles)
    fileCounter = fileCounter + 1;
    allFiles{fileCounter} = sbcFiles{ii};
    cellType = [cellType 4];
    matchFound = false;
    for jj = 1:ii-1
        if strcmp(sbcFiles{ii}(1:4), sbcFiles{jj}(1:4))
            matchFound = true;
            sameCells = [sameCells sameCells(jj+fileCountBtwnTypes)];
        end
    end
    if ~matchFound
        cellIDCounter = cellIDCounter + 1;
        sameCells = [sameCells cellIDCounter];
    end
end


%% loading and checking data

allDataTri = cell(size(allFiles));
allDataBi = cell(size(allFiles));
allCurvesTri = cell(size(allFiles));
allCurvesBi = cell(size(allFiles));
allThreshTri = zeros(size(allFiles));
allThreshBi = zeros(size(allFiles));
allThreshRatio = zeros(size(allFiles));

for i = 1:length(allFiles)
    % triphasic responses
    load(['data002/elecResp_' allFiles{i}])

    
    %checks to make elecResp analysis is current and correct
    elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', 0, 'keepLogBased', 0, 'plotInconsistentFits', true);

    allThreshTri(i) = elecResp.analysis.threshold; 

    save(['data002/elecResp_' allFiles{i}], 'elecResp')
    clear elecResp
    
    
    % biphasic responses
    load(['data003/elecResp_' allFiles{i}])

    %checks to make elecResp analysis is current and correct
    elecResp = checkForUnfinishedAnalysis(elecResp, 100, 'recalcAll', 0, 'keepLogBased', 0);

    allThreshBi(i) = elecResp.analysis.threshold; 
    
    allThreshRatio(i) = allThreshTri(i)/allThreshBi(i);
    
    save(['data003/elecResp_' allFiles{i}], 'elecResp')
    clear elecResp    
end


%% plotting

%colors = hsv(4);

colors = [204 0 0;
    50 70 247;
    90 156 0;
    214 149 68]/255;

markerTypes = {'x', '*', '+', 'o', 'd', 'h'};
%markerSizes = [12 12 12 8 8 8];
markerSizes = [8 8 8 6 6 6];
markerFill = [1 1 1 0 0 0];


%% "bar" graph


figure('position', [100 100 200 400])
hold on

currentMarkerIndices = sameCells;
switchIndex = 1;
for ii = 1:length(allFiles)
    if ii > 1 && cellType(ii) ~= cellType(ii-1) %switching between cell types
        
        meanThreshRatio = mean(allThreshRatio(switchIndex:ii-1));
        stdThreshRatio = std(allThreshRatio(switchIndex:ii-1));
        
        meanThreshCellType = cellType(ii-1);

        plot(0.2 + meanThreshCellType, meanThreshRatio, 's', 'MarkerEdgeColor', colors(meanThreshCellType,:), 'MarkerFaceColor', colors(meanThreshCellType,:))
        plot([0.2 + meanThreshCellType, 0.2 + meanThreshCellType], [meanThreshRatio - stdThreshRatio, meanThreshRatio + stdThreshRatio],...
            'color', colors(meanThreshCellType,:));
        
        switchIndex = ii;
        currentMarkerIndices = sameCells - (sameCells(ii)-1)*ones(size(sameCells));
        firstSameCell = sameCells(ii);
    end
    plotColor = colors(cellType(ii),:);
    markerType = markerTypes{currentMarkerIndices(ii)};
    markerSize = markerSizes(currentMarkerIndices(ii));
    
    if markerFill(currentMarkerIndices(ii))
        fillString = ['''MarkerFaceColor'', [' num2str(plotColor) '],'];
    else
    fillString = '';
    end
    
    plotCommand = ['plot(' num2str(cellType(ii)) ',' num2str(allThreshRatio(ii)) ',''' markerType,...
        ''',' '''MarkerEdgeColor'', [' num2str(plotColor) '] , ' fillString '''MarkerSize'', ' num2str(markerSize) ')'];
    eval(plotCommand)
end


set(gca, 'xlim', [0.5 4.5], 'ylim', [1 2], 'xtick', [])
text(2, 1.95, 'on parasol', 'color', colors(1,:))
text(2, 1.9, 'off parasol', 'color', colors(2,:))
text(2, 1.85, 'on midget', 'color', colors(3,:))
text(2, 1.8, 'SBC', 'color', colors(4,:))
ylabel('biphasic threshold / triphasic threshold')
%text(4, 1.8, 'axon', 'color', axonColor)


%% condensed bar graph
figure
hold on
for ii = 1:length(allFiles)
    plot(allThreshRatio(ii), 1, 'k.', 'MarkerSize', 12)
end

meanThreshRatio = mean(allThreshRatio);
stdThreshRatio = std(allThreshRatio);


plot(meanThreshRatio, 1.2, 'ks')
plot([meanThreshRatio - stdThreshRatio, meanThreshRatio + stdThreshRatio], [1.2, 1.2], 'k-')

hold off
set(gca, 'ylim', [0.5 1.5], 'xlim', [1 2], 'ytick', [])
title('biphasic threshold / triphasic threshold')









