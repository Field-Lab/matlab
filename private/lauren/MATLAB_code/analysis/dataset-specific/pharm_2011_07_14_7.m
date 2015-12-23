
clear all

cellIDsVis =   [233     437    709    711     408     601    556    631    397];
cellIDsStim =  [339     500    712    711     486     528    559    632    396];
stimElecs =    [17      30     50     43      28      41     38     43     27];
cellTypes =  {'onP', 'onP', 'onP', 'offP', 'offP', 'onM', 'onM', 'onM', 'offM'};
phaseWidth = 100;
pathToElecRespFilesBefore = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data002/';
pathToElecRespFilesAfter  = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data006/';

pathToNeuronsFileBefore = ['/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data004-005-no_refit/'...
    'data004-from-data004_data005-cf-lh/data004-from-data004_data005-cf-lh.neurons'];

pathToNeuronsFileAfter  = ['/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data004-005-no_refit/'...
    'data005-from-data004_data005-cf-lh/data005-from-data004_data005-cf-lh.neurons'];

rasterPadding = [0 0];

%% plot summary figure for each cell

for ii = 1:length(cellIDsStim)
    pathToElecResp.before = [pathToElecRespFilesBefore 'elecResp_n'...
        num2str(cellIDsStim(ii)) '_p' num2str(stimElecs(ii)) '_w' num2str(phaseWidth)];
    pathToElecResp.after  = [pathToElecRespFilesAfter  'elecResp_n'...
        num2str(cellIDsStim(ii)) '_p' num2str(stimElecs(ii)) '_w' num2str(phaseWidth)];
    
    pharm_results_plotter(cellIDsVis(ii), cellTypes{ii}, pathToElecResp,...
        pathToNeuronsFileBefore, 'pathToNeuronsFileAfter', pathToNeuronsFileAfter)
end




%% plot rasters for all cells

if 0
    for z = 1 %for code folding purposes only       
        
        
        % load visual responses %%%%%%%%%%%%%%%%%%%%%%%%%
        
        % before blockers added
        [spikeTimesAllWhite1 interval] = loadMovingBarResponses(pathToNeuronsFileBefore, cellIDsVis, 'white', 'trialPadding', rasterPadding);
        [spikeTimesAllBlack1]          = loadMovingBarResponses(pathToNeuronsFileBefore, cellIDsVis, 'black', 'trialPadding', rasterPadding);
        
        % after blockers added
        [spikeTimesAllWhite2]          = loadMovingBarResponses(pathToNeuronsFileAfter, cellIDsVis, 'white', 'trialPadding', rasterPadding);
        [spikeTimesAllBlack2]          = loadMovingBarResponses(pathToNeuronsFileAfter, cellIDsVis, 'black', 'trialPadding', rasterPadding);
        
        nCells = length(cellIDsVis);
        
        nBarReps1 = length(spikeTimesAllWhite1);
        nBarReps2 = length(spikeTimesAllWhite2);
        
        lInterval = interval{1}(2) - interval{1}(1);
        
        gapSize = 10;
        
        
        figure
        for i = 1:nCells
            axes('position', [0.1 (nCells - i + 1)/(nCells+1.5) 0.35 1/(nCells+2.5)])
            hold on
            fill([0 lInterval lInterval 0], [nBarReps1 nBarReps1 nBarReps1+gapSize nBarReps1+gapSize],...
                [0.6 0.6 0.6], 'EdgeColor', 'none')
            fill([0 lInterval lInterval 0], [nBarReps1+gapSize nBarReps1+gapSize nBarReps1+gapSize+nBarReps2 nBarReps1+gapSize+nBarReps2],...
                [0.95 0.95 0.95], 'EdgeColor', 'none')
            for k = 1:nBarReps1
                for j = 1:length(spikeTimesAllWhite1{k}{i})
                    plot([spikeTimesAllWhite1{k}{i}(j) spikeTimesAllWhite1{k}{i}(j)], [k-1 k], 'k-', 'LineWidth', 1)
                end
            end
            
            for k = 1:nBarReps2
                for j = 1:length(spikeTimesAllWhite2{k}{i})
                    plot([spikeTimesAllWhite2{k}{i}(j) spikeTimesAllWhite2{k}{i}(j)], [k-1 k]+nBarReps1+gapSize, 'k-', 'LineWidth', 1)
                end
            end

            hold off
            
            set(gca, 'yLim', [0 nBarReps1+nBarReps2+gapSize], 'xlim', [0 lInterval], 'yDir', 'reverse')
            if i == 1
                title('white bar')
            end
            if i == nCells
                xlabel('time (s)')
            else
                set(gca, 'xtick', [])
            end
            ylabel(['cell' num2str(cellIDsStim(i))])
            
            axes('position', [0.55 (nCells - i + 1)/(nCells+1.5) 0.35 1/(nCells+2.5)])
            hold on
            fill([0 lInterval lInterval 0], [nBarReps1 nBarReps1 nBarReps1+gapSize nBarReps1+gapSize],...
                [0.6 0.6 0.6], 'EdgeColor', 'none')
            fill([0 lInterval lInterval 0], [nBarReps1+gapSize nBarReps1+gapSize nBarReps1+gapSize+nBarReps2 nBarReps1+gapSize+nBarReps2],...
                [0.95 0.95 0.95], 'EdgeColor', 'none')
            for k = 1:nBarReps1
                for j = 1:length(spikeTimesAllBlack1{k}{i})
                    plot([spikeTimesAllBlack1{k}{i}(j) spikeTimesAllBlack1{k}{i}(j)], [k-1 k], 'k-', 'LineWidth', 1)
                end
            end
            for k = 1:nBarReps2
                for j = 1:length(spikeTimesAllBlack2{k}{i})
                    plot([spikeTimesAllBlack2{k}{i}(j) spikeTimesAllBlack2{k}{i}(j)], [k-1 k]+nBarReps1+gapSize, 'k-', 'LineWidth', 1)
                end
            end
            
            hold off
            
            set(gca, 'yLim', [0 nBarReps1+nBarReps2+gapSize], 'xlim', [0 lInterval], 'yDir', 'reverse')
            if i == 1
                title('black bar')
            end
            if i == nCells
                xlabel('time (s)')
            else
                set(gca, 'xtick', [])
            end
        end
    end
end
