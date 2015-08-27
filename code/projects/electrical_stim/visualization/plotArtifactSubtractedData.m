function plotArtifactSubtractedData(pathToElecRespFile,nameOfElecRespFile, amplitudeRange)
% PLOTARTIFACTSUBTRACTEDDATA(...) plots
% Inputs:       pathToElecRespFile: string
%               nameOfElecRespFile: string
%               amplitudeRange: vector [size(amplitudeRange) = (1,2)] specifying lower and upper bounds on the stimulation amplitudes that will be plotted       
% Usage: plotArtifactSubtractedData('/Volumes/Analysis/2015-04-14-0/data001/','elecResp_n3947_p264.mat', [0.6 1] )
% Lauren Grosberg 8/2015

figure; set(gcf,'Position', [1         191        1907         907]);

% load elecResp files
filename = fullfile(pathToElecRespFile,nameOfElecRespFile); 
try
    temp = load(filename);
catch
    disp([filename ' not found']);
    return;
end
elecResp = temp.elecResp; clear temp; 

% Get response curves
nMovies = length(elecResp.stimInfo.movieNos);
data = zeros(2, nMovies);
data(1,:) = elecResp.stimInfo.stimAmps;
data(2,:) = elecResp.analysis.successRates;
data(3,:) = elecResp.stimInfo.nPulses;
lockedAmps = elecResp.analysis.finalized;
for i = length(elecResp.stimInfo.stimAmps): -1: 1
    if isempty(elecResp.analysis.type{i})
        data(:,i) = [];
        lockedAmps(i) = [];
    end
end

% linear-based
data(1,:) = abs(data(1,:));
[erfParams, completeFit, ~] = erfFitter(data, 2, -1, 'makePlot', 0, 'lockedAmps', lockedAmps);
thresholdHum = -erfParams(2)/erfParams(1);
curveHum = completeFit; 




%% Plot the raw data traces - artifact estimated by the human
centerChannel = elecResp.cells.recElec;
firstMovieIdx = find(abs(elecResp.stimInfo.stimAmps) > amplitudeRange(1),1,'first');
lastMovieIdx =  find(abs(elecResp.stimInfo.stimAmps) < amplitudeRange(2),1,'last');
moviesToPlot = firstMovieIdx:lastMovieIdx; 
numSubplots = length(moviesToPlot) + 1; 
for ii = 1:length(moviesToPlot)
    movieIndex = moviesToPlot(ii);
    if length(unique(elecResp.stimInfo.stimAmps)) == length(elecResp.stimInfo.stimAmps) % There are not 'movies' with duplicate stimulation amps
        subplot(1, numSubplots, 1+ii);
        dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
            elecResp.stimInfo.movieNos(movieIndex), 99999);
        subtractionVector = elecResp.analysis.estArtifact{movieIndex}(elecResp.cells.goodElecs == centerChannel, :)';
        dataToPlot = zeros(size(dataTraces, 1), size(dataTraces, 3));
        nPulses = size(dataTraces, 1);
        for i = 1:nPulses
            dataToPlot(i, :) = squeeze(dataTraces(i, centerChannel, :)) - subtractionVector;
        end
        latencies = [elecResp.analysis.latencies{movieIndex} elecResp.analysis.otherLatencies{movieIndex}];
        
        % Plot neuron template.
        neuronTemplate = elecResp.cells.mainEI(centerChannel,:);
        alignmentPoint = find(neuronTemplate == min(neuronTemplate));
        hold on;
        for i = 1:elecResp.stimInfo.nPulses(movieIndex)
            if ~any(latencies(i,:))
                misses = plot(dataToPlot(i, :), 'Color',[0.8 0.8 0.8]);
            end
        end
        for i = 1:elecResp.stimInfo.nPulses(movieIndex)
            if any(latencies(i,:))
                shiftedData = circshift(dataToPlot(i, :)', alignmentPoint - round(latencies(i,:)));
                hits = plot(shiftedData, 'r');
            end
        end
        plot(neuronTemplate,'k','LineWidth',1);
        ylim([min(neuronTemplate)-50 max(neuronTemplate)+50])
        xlim([0 40]);
        title(sprintf('human\n%0.0f: stimAmp %0.3f uA', movieIndex,...
            elecResp.stimInfo.stimAmps(movieIndex)));
    else % There are movies with duplicate stimulation amps.
        keyboard; % figure out how to handle this still. 
%         subplot(2, 5, 1+ii);
%         % Read in all the info from the particular stimulation amp
%         if ii == 1
%             adjustedIndex = find(elecResp.stimInfo.stimAmps == elecResp.stimInfo.stimAmps(movieIndex));
%         else
%             adjustedIndex = adjustedIndex + 2;
%         end
%         
%         for aa = 1:length(adjustedIndex)
%             dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
%                 elecResp.stimInfo.movieNos(adjustedIndex(aa)), 99999);
%             subtractionVector = elecResp.analysis.estArtifact{adjustedIndex(aa)}(elecResp.cells.goodElecs == centerChannel, :)';
%             latencies = [elecResp.analysis.latencies{adjustedIndex(aa)} elecResp.analysis.otherLatencies{adjustedIndex(aa)}];
%             dataToPlot = zeros(size(dataTraces, 1), size(dataTraces, 3));
%             nPulses = size(dataTraces, 1);
%             for i = 1:nPulses
%                 dataToPlot(i, :) = squeeze(dataTraces(i, centerChannel, :)) - subtractionVector;
%             end
%             
%             % Plot neuron template.
%             neuronTemplate = elecRespAuto.neuronInfo.templates{1};
%             alignmentPoint = find(neuronTemplate == min(neuronTemplate));
%             hold on;
%             for i = 1:elecResp.stimInfo.nPulses(movieIndex)
%                 if ~any(latencies(i,:))
%                     misses = plot(dataToPlot(i, :), 'Color',[0.8 0.8 0.8]);
%                 else
%                     misses = plot(1, 'Color',[0.8 0.8 0.8]);
%                 end
%             end
%             for i = 1:elecResp.stimInfo.nPulses(movieIndex)
%                 if any(latencies(i,:))
%                     shiftedData = circshift(dataToPlot(i, :)', alignmentPoint - round(latencies(i,:)));
%                     hits = plot(shiftedData, 'r');
%                 else
%                     hits = plot(1,'r');
%                 end
%             end
%         end
%         
%         plot(neuronTemplate,'k','LineWidth',1);
%         ylim([min(neuronTemplate)-50 max(neuronTemplate)+50])
%         xlim([0 40]);
%         title(sprintf('human\n%0.0f&%0.0f: stimAmp %0.3f uA', adjustedIndex(1),...
%             adjustedIndex(2),elecResp.stimInfo.stimAmps(adjustedIndex(1))));
    end
    
end
try legend([hits misses],'spikes', 'no spikes'); 
catch
end

% Plot the erf fits
subplot(1,numSubplots,1); 
colors = lines(2);
plot(abs(elecResp.stimInfo.stimAmps),elecResp.analysis.successRates,...
    'o','LineWidth',3,'Color',colors(1,:)); grid on; 
hold on; plot(curveHum(1,:),curveHum(2,:),'LineWidth',3,'Color',colors(1,:)); 
title(sprintf('erf fits human thresh: %0.02f',thresholdHum)); 