%% This script plots single electrode stimulation of a rat retina during TTX perfusion
% The purpose is to observe the variability in the artifact between
% electrodes and for different stimulation amplitudes for 1-electrode stim
% Dataset 2014-07-24-0
%% Plot the artifact recorded using TTX
colors = colormap(jet);
colors = colors(1:2:end,:);
saveFiles = 0;
for p = [40 36 30 29 9 380 371 363 374 325 381 347];
    pathToAnalysisData = '/Volumes/Analysis/2014-07-24-0/data017/';
    patternNo          = p;
    movieNos           = findMovieNos(pathToAnalysisData,patternNo);
    % recordingElec      = elecResp.cells.recElec;
    recordingElec      = patternNo;
    sampRate           = 20000;
    textPositions = sort(linspace(-1200,400,size(movieNos,2)));
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
          movieNos(1), 99999);
    lastTraces = mean(squeeze(dataTraces(:,recordingElec,:)));
    figure;
    for m = 2:size(movieNos, 2)
        
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
            movieNos(m), 99999);
        theseTraces = mean(squeeze(dataTraces(:,recordingElec,:)));

        tracesToPlot = theseTraces; 
        timepoints = linspace(0,size(tracesToPlot,2)/sampRate * 1000, size(tracesToPlot,2));
        hold on; plot(timepoints, tracesToPlot','Color',colors(m,:));
        
        stimAmp = getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
        text(5.1,textPositions(m),num2str(stimAmp),'Color',colors(m,:));
        xlim([0 6]); ylim([-1200 500]);
        lastTraces = theseTraces; 
    end
    
    xlabel('Time (ms)');
    title(sprintf('%s pattern %d \nTTX recordings',pathToAnalysisData,patternNo));
    
    if saveFiles
        savename= ['/Users/grosberg/Desktop/forGonzalo/TTX_traces_201407240_data017_p' num2str(patternNo)];
        %     saveas(gcf,'test','fig');
        hgsave(gcf,[savename '.fig']);
        %      saveas(gcf,eval(savename));
        % savefig(savename,'eps');
    end
end