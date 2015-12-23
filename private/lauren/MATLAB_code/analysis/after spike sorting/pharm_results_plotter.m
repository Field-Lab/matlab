function pharm_results_plotter(cellIDVis, cellType, pathToElecResp, pathToNeuronsFile, varargin)


p = inputParser;

p.addRequired('cellIDVis', @isnumeric)
p.addRequired('cellType', @ischar)
p.addRequired('pathToElecResp', @isstruct) %should have fields .before and .after
p.addRequired('pathToNeuronsFile', @ischar)


p.addParamValue('rasterPadding', [0 0], @isnumeric)
p.addParamValue('pathToNeuronsFileAfter', [], @ischar) %use if moving bar datasets of before/after drugs are separate

p.parse(cellIDVis, cellType, pathToElecResp, pathToNeuronsFile, varargin{:})

rasterPadding = p.Results.rasterPadding;
pathToNeuronsFileAfter = p.Results.pathToNeuronsFileAfter;


fileseps = strfind(pathToElecResp.before, '/');
elecRespName = pathToElecResp.before(fileseps(end)+1:end);


%% load elecrical responses

elecResp1 = load(pathToElecResp.before);
elecResp2 = load(pathToElecResp.after);

elecResp1 = checkForUnfinishedAnalysis(elecResp1.elecResp, 100);
elecResp2 = checkForUnfinishedAnalysis(elecResp2.elecResp, 100);

erfParams1 = elecResp1.analysis.erfParams;
erfParams2 = elecResp2.analysis.erfParams;

finalized1 = elecResp1.analysis.finalized;
finalized2 = elecResp2.analysis.finalized;

succRates1 = elecResp1.analysis.successRates;
succRates2 = elecResp2.analysis.successRates;

stimAmps1 = elecResp1.stimInfo.stimAmps;
stimAmps2 = elecResp2.stimInfo.stimAmps;

clear elecResp1 elecResp2

%% load visual responses

if isempty(pathToNeuronsFileAfter) %same dataset for before and after drug application
    [spikeTimesAllWhite interval] = loadMovingBarResponses(pathToNeuronsFile, cellIDVis, 'white', 'trialPadding', rasterPadding);
    [spikeTimesAllBlack]          = loadMovingBarResponses(pathToNeuronsFile, cellIDVis, 'black', 'trialPadding', rasterPadding);
    
    nBarReps = length(spikeTimesAllWhite);

else %separate moving bar datasets for before and after drug application
    % before blockers added
    [spikeTimesAllWhite1 interval] = loadMovingBarResponses(pathToNeuronsFile, cellIDVis, 'white', 'trialPadding', rasterPadding);
    [spikeTimesAllBlack1]          = loadMovingBarResponses(pathToNeuronsFile, cellIDVis, 'black', 'trialPadding', rasterPadding);
    
    % after blockers added
    [spikeTimesAllWhite2]          = loadMovingBarResponses(pathToNeuronsFileAfter, cellIDVis, 'white', 'trialPadding', rasterPadding);
    [spikeTimesAllBlack2]          = loadMovingBarResponses(pathToNeuronsFileAfter, cellIDVis, 'black', 'trialPadding', rasterPadding);
        
    nBarReps1 = length(spikeTimesAllWhite1);
    nBarReps2 = length(spikeTimesAllWhite2);
    
    gapSize = 10;
end

lInterval = interval{1}(2) - interval{1}(1);

%% summary plot

figure('position', [100 100 800 400])

%response curves
axes('position', [0.55 0.1 0.35 0.8])
hold on

xProj1 = min(abs(stimAmps1(~~finalized1))):0.001:max(abs(stimAmps1(~~finalized1)));
xProj2 = min(abs(stimAmps2(~~finalized2))):0.001:max(abs(stimAmps2(~~finalized2)));

projection1 = 0.5 + 0.5*erf(erfParams1(1)*xProj1+erfParams1(2));
projection2 = 0.5 + 0.5*erf(erfParams2(1)*xProj2+erfParams2(2));

plot(xProj1, projection1, 'r')
plot(xProj2, projection2, 'b')

for jj = 1:length(stimAmps1)
    if finalized1(jj)
        plot(abs(stimAmps1(jj)), succRates1(jj), 'ro')
    end
end
for jj = 1:length(stimAmps2)
    if finalized2(jj)
        plot(abs(stimAmps2(jj)), succRates2(jj), 'bo')
    end
end

hold off
title([elecRespName ' (' cellType ')'])

% rasters
if isempty(pathToNeuronsFileAfter) %same dataset for before and after drug application

    axes('position', [0.1 0.1 0.35 0.35])
    hold on
    
%     fill([0 lInterval lInterval 0], [nBarReps1 nBarReps1 nBarReps1+gapSize nBarReps1+gapSize],...
%         [0.6 0.6 0.6], 'EdgeColor', 'none')
%     fill([0 lInterval lInterval 0], [nBarReps1+gapSize nBarReps1+gapSize nBarReps1+gapSize+nBarReps2 nBarReps1+gapSize+nBarReps2],...
%         [0.95 0.95 0.95], 'EdgeColor', 'none')
    for k = 1:nBarReps
        for j = 1:length(spikeTimesAllWhite{k}{1})
            plot([spikeTimesAllWhite{k}{1}(j) spikeTimesAllWhite{k}{1}(j)], [k-1 k], 'k-', 'LineWidth', 1)
        end
    end
    
    hold off
    set(gca, 'yLim', [0 nBarReps], 'xlim', [0 lInterval], 'yDir', 'reverse')
    title('white bar')
    xlabel('time (s)')
    
    axes('position', [0.1 0.55 0.35 0.35])
    hold on
%     fill([0 lInterval lInterval 0], [nBarReps1 nBarReps1 nBarReps1+gapSize nBarReps1+gapSize],...
%         [0.6 0.6 0.6], 'EdgeColor', 'none')
%     fill([0 lInterval lInterval 0], [nBarReps1+gapSize nBarReps1+gapSize nBarReps1+gapSize+nBarReps2 nBarReps1+gapSize+nBarReps2],...
%         [0.95 0.95 0.95], 'EdgeColor', 'none')
    for k = 1:nBarReps
        for j = 1:length(spikeTimesAllBlack{k}{1})
            plot([spikeTimesAllBlack{k}{1}(j) spikeTimesAllBlack{k}{1}(j)], [k-1 k], 'k-', 'LineWidth', 1)
        end
    end
    
    hold off
    set(gca, 'yLim', [0 nBarReps], 'xlim', [0 lInterval], 'yDir', 'reverse')
    title('black bar')
    
else %separate moving bar datasets for before and after drug application

    axes('position', [0.1 0.1 0.35 0.35])
    hold on
    
    fill([0 lInterval lInterval 0], [nBarReps1 nBarReps1 nBarReps1+gapSize nBarReps1+gapSize],...
        [0.6 0.6 0.6], 'EdgeColor', 'none')
    fill([0 lInterval lInterval 0], [nBarReps1+gapSize nBarReps1+gapSize nBarReps1+gapSize+nBarReps2 nBarReps1+gapSize+nBarReps2],...
        [0.95 0.95 0.95], 'EdgeColor', 'none')
    for k = 1:nBarReps1
        for j = 1:length(spikeTimesAllWhite1{k}{1})
            plot([spikeTimesAllWhite1{k}{1}(j) spikeTimesAllWhite1{k}{1}(j)], [k-1 k], 'k-', 'LineWidth', 1)
        end
    end
    
    for k = 1:nBarReps2
        for j = 1:length(spikeTimesAllWhite2{k}{1})
            plot([spikeTimesAllWhite2{k}{1}(j) spikeTimesAllWhite2{k}{1}(j)], [k-1 k]+nBarReps1+gapSize, 'k-', 'LineWidth', 1)
        end
    end
    
    hold off
    set(gca, 'yLim', [0 nBarReps1+nBarReps2+gapSize], 'xlim', [0 lInterval], 'yDir', 'reverse')
    title('white bar')
    xlabel('time (s)')
    
    axes('position', [0.1 0.55 0.35 0.35])
    hold on
    fill([0 lInterval lInterval 0], [nBarReps1 nBarReps1 nBarReps1+gapSize nBarReps1+gapSize],...
        [0.6 0.6 0.6], 'EdgeColor', 'none')
    fill([0 lInterval lInterval 0], [nBarReps1+gapSize nBarReps1+gapSize nBarReps1+gapSize+nBarReps2 nBarReps1+gapSize+nBarReps2],...
        [0.95 0.95 0.95], 'EdgeColor', 'none')
    for k = 1:nBarReps1
        for j = 1:length(spikeTimesAllBlack1{k}{1})
            plot([spikeTimesAllBlack1{k}{1}(j) spikeTimesAllBlack1{k}{1}(j)], [k-1 k], 'k-', 'LineWidth', 1)
        end
    end
    
    for k = 1:nBarReps2
        for j = 1:length(spikeTimesAllBlack2{k}{1})
            plot([spikeTimesAllBlack2{k}{1}(j) spikeTimesAllBlack2{k}{1}(j)], [k-1 k]+nBarReps1+gapSize, 'k-', 'LineWidth', 1)
        end
    end
    hold off
    set(gca, 'yLim', [0 nBarReps1+nBarReps2+gapSize], 'xlim', [0 lInterval], 'yDir', 'reverse')
    title('black bar')

end

