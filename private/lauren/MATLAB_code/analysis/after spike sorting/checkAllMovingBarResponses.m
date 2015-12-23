%% load visual responses

clear all

%pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2011-08-01-0/data011/data011.neurons'; %path to file with drugs applied
pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2011-08-01-0/data011-013-no_refit/data011-013-no_refit.neurons';

datarun.names.rrs_neurons_path = pathToNeuronsFile;
datarun = load_neurons(datarun);

[spikeTimesAllWhite interval] = loadMovingBarResponses(pathToNeuronsFile, datarun.cell_ids, 'white');
[spikeTimesAllBlack]          = loadMovingBarResponses(pathToNeuronsFile, datarun.cell_ids, 'black');
    
nBarReps = min(length(spikeTimesAllWhite), length(spikeTimesAllBlack));
lInterval = interval{1}(2) - interval{1}(1);



%% plot all visual responses in a slider plot

% makes slider plot: needs to be updated to reflect log vs. not-log
sliderFig = figure;
slider = make_loop_slider_list(1,1,length(datarun.cell_ids));

topAxes = axes('position', [0.1 0.55 0.8 0.35]);
bottomAxes = axes('position', [0.1 0.1 0.8 0.35]);
titleAxes = axes('position', [0.1 0.95 0.8 0.05]);
axis off


while ishandle(sliderFig)
    ii = round(get(slider,'Value'));
    %cla
    
    cla(titleAxes)
    axes(titleAxes)
    text(0,0, ['cell ' num2str(datarun.cell_ids(ii))])
    axis off
    
    cla(topAxes)
    cla(bottomAxes)
    
    drawnow

    axes(topAxes)
    hold on
    for k = 1:nBarReps
        for j = 1:length(spikeTimesAllBlack{k}{ii})
            plot([spikeTimesAllBlack{k}{ii}(j) spikeTimesAllBlack{k}{ii}(j)], [k-1 k], 'k-', 'LineWidth', 1)
        end
    end
    hold off
    set(gca, 'yLim', [0 nBarReps], 'xlim', [0 lInterval], 'yDir', 'reverse')
    title('black bar')
    
    drawnow
    
    axes(bottomAxes)
    hold on
    for k = 1:nBarReps
        for j = 1:length(spikeTimesAllWhite{k}{ii})
            plot([spikeTimesAllWhite{k}{ii}(j) spikeTimesAllWhite{k}{ii}(j)], [k-1 k], 'k-', 'LineWidth', 1)
        end
    end
    hold off
    set(gca, 'yLim', [0 nBarReps], 'xlim', [0 lInterval], 'yDir', 'reverse')
    title('white bar')
    xlabel('time (s)')

    uiwait;
end
