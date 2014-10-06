% Analyze pulse triplet rat experiment from 2014-07-24-0 to determine if a
% large pre-pulse can silence all neural activity in a region and therefore
% be used for artifact estimation

%% Plot the artifact & activation curves for a single cell and pattern no.

pathToAnalysisData = '/Volumes/Analysis/2014-07-24-0/data003/';
patternNo          = 347;
neuron             = 5267;
fname = [pathToAnalysisData 'elecResp_n' num2str(neuron) '_p' num2str(patternNo) '.mat'];

if exist(fname,'file') == 2
    tmp       = load(fname);
    elecResp = tmp.elecResp;
else
    errordlg([fname 'does not exist']);
end
colors = colormap(lines);
colors = colors(1:2:end,:);
figure;
for ii = 1:1:15 %size(elecResp.analysis.estArtifact,1)
    art = elecResp.analysis.estArtifact{ii};
    hold on; plot(art,'Color',colors(ii,:));
end
title(sprintf('Artifact estimation using analysisGUI spike sorting \nstim amps from %0.2f to %0.2f uA \nrecording elec %d, neuron %d',elecResp.stimInfo.stimAmps(1),elecResp.stimInfo.stimAmps(ii),elecResp.cells.recElec,neuron));
%% Plot the artifact recorded using TTX
% for p = [40 36 30 29 9 380 371 363 374 325 381];
    pathToAnalysisData = '/Volumes/Analysis/2014-07-24-0/data017/';
    patternNo          = 347; %p;
    movieNos           = findMovieNos(pathToAnalysisData,patternNo);
    % recordingElec      = elecResp.cells.recElec;
    recordingElec      = patternNo;
    sampRate           = 20000;
    textPositions = sort(linspace(-1200,400,size(movieNos,2)));
    
    figure;
    for m = 1:size(movieNos, 2)
        
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
            movieNos(m), 99999);
        
        tracesToPlot = squeeze(dataTraces(:,recordingElec,:)); %trials, electrodes, timepoints
        timepoints = linspace(0,size(tracesToPlot,2)/sampRate * 1000, size(tracesToPlot,2));
        hold on; plot(timepoints, tracesToPlot','Color',colors(m,:));
        
        stimAmp = getStimAmps(pathToAnalysisData, patternNo, movieNos(m), 512);
        text(5.1,textPositions(m),num2str(stimAmp),'Color',colors(m,:));
        xlim([0 6]); ylim([-1200 500]);
    end
    
    xlabel('Time (ms)');
    title(sprintf('%s pattern %d \nTTX recordings',pathToAnalysisData,patternNo));
    
    savename= ['/Users/grosberg/Desktop/temp files/TTX_traces_2014-07-24-0-data017-p' num2str(patternNo) ];
%     saveas(gcf,savename,'jpg');
        
% end
        % For cell number 5267, plot the estimated artifact.
        
        % filePath  = '/Volumes/Analysis/2014-07-24-0/data006/pattern_files/';
        % tic
        % elecsUsed = getElecsFromPatternFiles(filePath);
        % toc
        % To look at the increased amplitude, add 512 to the electrode number
        %% Plot the artifact using the pulse triplets
        
        pathToAnalysisData = '/Volumes/Analysis/2014-07-24-0/data006/'; % Pulse triplet experiment, low stimulation amps.
        elecNum            = 347;
        patternNo          = elecNum+512*4;
        movieNo            = 77;
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNo, 99999);
        
        figure; plot(timepoints, tracesToPlot');
        xlabel('Time (ms)');
        title(sprintf('%s pattern %d; movie # %d',pathToAnalysisData,patternNo,movieNo));
        
        % Find movie nos.
        movieNos      = findMovieNos(pathToAnalysisData,patternNo);
        %%
        pathToAnalysisData = '/Volumes/Analysis/2014-07-24-0/long-chunks/data006';
        % pathToAnalysisData = '/Volumes/Analysis/2014-07-24-0/data006/'; % Pulse triplet experiment, low stimulation amps.
        % patterns           = [347 859 1371 1883 2395 2907 3419];
        
        sampRate          = 20000;
        
        patternNo         = 225; %[347 859 1371];
        recordingElec     = patternNo;
        movieNos          = findMovieNos(pathToAnalysisData,patternNo);
        
        figure;
        for m = 1:length(movieNos)
        movieNo      = movieNos(m);
        dataTraces   = NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
            movieNo, 99999);
        tracesToPlot = squeeze(dataTraces(:,recordingElec,:)); %trials, electrodes, timepoints
        timepoints = linspace(0,size(tracesToPlot,2)/sampRate * 1000, size(tracesToPlot,2));
        plot(timepoints, tracesToPlot');
        title(sprintf('%s \npattern %d; movie # %d',pathToAnalysisData,patternNo,movieNo));
        xlabel('Time (ms)');
        pause;
        end
        %%
        figure;
        for m = 1:length(movieNos)
            movieNo = movieNos(m);
            dataTraces   = NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
                movieNo, 99999);
            tracesToPlot = squeeze(dataTraces(:,recordingElec,:)); %trials, electrodes, timepoints
            subplot(2,3,m);
            plot(timepoints, tracesToPlot');
            xlabel('Time (ms)');
            title(sprintf('%s pattern %d; movie # %d',pathToAnalysisData,patternNo,movieNo));
        end
        movieNo       = movieNos(1);
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
            movieNo, 99999);
        tracesToPlot = squeeze(dataTraces(:,recordingElec,:)); %trials, electrodes, timepoints
        largeAmp1 = 347;
        largeAmp2 = 859;
        figure;
        for p = 1:length(patterns)
            patternNo     = patterns(p);
            movieNos      = findMovieNos(pathToAnalysisData,patternNo);
            %     disp(movieNos);
            movieNo       = movieNos(1);
            dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
                movieNo, 99999);
            tracesToPlot = squeeze(dataTraces(:,recordingElec,:)); %trials, electrodes, timepoints
            timepoints = linspace(0,size(tracesToPlot,2)/sampRate * 1000, size(tracesToPlot,2));
            
            keyboard;
            plot(timepoints, tracesToPlot');
            xlabel('Time (ms)');
            title(sprintf('%s pattern %d; movie # %d',pathToAnalysisData,patternNo,movieNo));
            pause;
        end
        
