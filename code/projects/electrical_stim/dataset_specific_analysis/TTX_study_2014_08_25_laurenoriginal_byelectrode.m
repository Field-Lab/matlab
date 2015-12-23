%% This script plots single electrode stimulation of a rat retina during TTX perfusion
% The purpose is to observe the variability in the artifact between
% electrodes and for different stimulation amplitudes for 1-electrode stim
% Dataset 2014-07-24-0
%% Plot the artifact recorded using TTX

%Initial params
pathToAnalysisData = '/Volumes/Analysis/2014-07-24-0/data017/';
saveFigs = 0;

%get all movieNos
pcnt = 1; 
movieNosCell = {}; 
patts=[40 36 30 29 9 380 371 363 374 325 381 347]; 
for p = patts;
	movieNosCell{pcnt} =  findMovieNos(pathToAnalysisData,p); 
	pcnt = pcnt + 1;
end

%get minimum # of movies for any one pattern
minMovNum = 1000000;
for i=1:pcnt-1;
	len = length(movieNosCell{i});
	if len < minMovNum; minMovNum = len; end
end
	
%define colormap
colors = colormap(jet);
colors = colors(1:2:end,:);
colors=colors(round(linspace(1,size(colors,1),length(patts))),:)

%Loop through movies
for m = 1:minMovNum;

    figure;
pcnt = 1;
for patternNo = [40 36 30 29 9 380 371 363 374 325 381 347];
    movieNo           = movieNosCell{pcnt}(m);
    % recordingElec      = elecResp.cells.recElec;
    recordingElec      = patternNo;
    sampRate           = 20000;
    textPositions = sort(linspace(-1200,400,size(movieNosCell{pcnt},2)));
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
          movieNo, 99999);
	theseTraces = mean(squeeze(dataTraces(:,recordingElec,:)));

        tracesToPlot = theseTraces; 
        timepoints = linspace(0,size(tracesToPlot,2)/sampRate * 1000, size(tracesToPlot,2));
        plot(timepoints, tracesToPlot','Color',colors(pcnt,:)); hold on; 
        
        stimAmp = getStimAmps(pathToAnalysisData, patternNo, movieNo);
        text(5.1,textPositions(m),num2str(stimAmp),'Color',colors(pcnt,:));
        xlim([0 6]); ylim([-1200 500]);
        lastTraces = theseTraces; 
    
    xlabel('Time (ms)');
    title(sprintf('%s pattern %d \nTTX recordings',pathToAnalysisData,patternNo));
	pcnt = pcnt + 1;
end
    
    if saveFiles
        savename= ['/Users/grosberg/Desktop/forGonzalo/TTX_traces_201407240_data017_p' num2str(patternNo)];
        %     saveas(gcf,'test','fig');
        hgsave(gcf,[savename '.fig']);
        %      saveas(gcf,eval(savename));
        % savefig(savename,'eps');
    end
end
