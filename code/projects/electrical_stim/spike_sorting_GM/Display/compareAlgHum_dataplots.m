function compareAlgHum_dataplots(elecRespAuto, n)
% Compares activation curves calculated by human sorting with those
% calculated using Gonzalo's spike sorting algorithm 
% there must exist an elecResp file in order to run this function
% Inputs:    elecRespAuto
%            n: index of the neuronId that should be plotted. n=1 if
%            there is only one neuron specified to GM algorithm

% Lauren Grosberg 8/2015

figure; set(gcf,'Position', [1         191        1907         907]);
neuronIds = elecRespAuto.neuronInfo.neuronIds;
neuronId = neuronIds(n);
e = elecRespAuto.neuronInfo.prefElectrodes{n}(1);
breakAxon    =  elecRespAuto.tracesInfo.breakAxon{e};
breakRecElec = elecRespAuto.tracesInfo.breakRecElecs{e};
J = elecRespAuto.tracesInfo.J;
I = elecRespAuto.tracesInfo.I;

spikeProbs    = nansum(elecRespAuto.spikes{n}')./I;
spikeLogProbs = elecRespAuto.LogisticReg(n,:);
amps = abs(elecRespAuto.stimInfo.listAmps(:,1))';


clear latencies
latencies = zeros(J,3);
for j=1:J
    lats = elecRespAuto.latencies{n}(j,1:I(j));
    lats = lats(lats>0);
    if(isempty(lats))
        latencies(j,1:3)=NaN;
    else
        
        latencies(j,1) = prctile(lats,25);
        latencies(j,2) = prctile(lats,50);
        latencies(j,3) = prctile(lats,75);
    end
end

latencies = latencies/elecRespAuto.params.sampRate*1000;

if ~strcmp(elecRespAuto.path(1:18),'/Volumes/Analysis/')
    gm_path = elecRespAuto.path;
    ii = find(gm_path==filesep,3,'last');
    my_path = fullfile('/Volumes/Analysis/',gm_path(ii(1):ii(3)));
    elecRespAuto.path = my_path;
end

fname = ['elecResp_n' num2str(neuronId) '_p' num2str(elecRespAuto.stimInfo.patternNo) '.mat']; 
filename = fullfile(elecRespAuto.path,fname); 
temp = load(filename); 
elecResp = temp.elecResp; 

subplot(2,5,1); 
hum = plot(abs(elecResp.stimInfo.stimAmps),elecResp.analysis.successRates,...
    'o-','LineWidth',3); grid on; 
hold on; alg = plot(amps,spikeProbs,'^-','linewidth',3); 
logRegAlg = plot(amps,spikeLogProbs,'--','linewidth',2); 
legend([hum alg(1,1) logRegAlg],'Human','Algorithm',...
    'Logistic Regression (algorithm)'); 
xlabel('stimulation amplitude (\muA)'); 
if length(elecResp.stimInfo.electrodes) == 2
    title(sprintf('n%0.0f p%0.0f\nstimE: %0.0f %0.0f recE: %0.0f\ntemplate: %s\nstim data: ..%s',...
        neuronId, elecResp.stimInfo.patternNo, elecResp.stimInfo.electrodes,elecResp.cells.recElec,...
        elecResp.names.rrs_ei_path(end-9:end),elecResp.names.data_path));
else
    title(sprintf('n%0.0f p%0.0f\nstimE: %0.0f recE: %0.0f\ntemplate: %s\nstim data: ..%s',...
        neuronId,elecResp.stimInfo.patternNo, elecResp.stimInfo.electrodes,elecResp.cells.recElec,...
        elecResp.names.rrs_ei_path(end-9:end),elecResp.names.data_path(end-21:end)));
end
ylim([0 1]); 
xlabel('stimulation amplitude (\muA)'); 
ylabel('spike probability'); 

for b = breakAxon
    plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','yellow')
end
axis([amps(1) amps(end) 0 1]); 


subplot(2,5,6); 
colors = lines(2);
[thresholdHum, thresholdAlg, curveHum, curveAlg] = fitToErfOutputAndHuman(elecRespAuto) ; 
plot(abs(elecResp.stimInfo.stimAmps),elecResp.analysis.successRates,...
    'o','LineWidth',3,'Color',colors(1,:)); grid on; 
hold on; plot(curveHum(1,:),curveHum(2,:),'LineWidth',3,'Color',colors(1,:)); 
hold on; plot(amps,spikeProbs,'^','linewidth',3,'Color', colors(2,:)); 
hold on; plot(curveAlg(1,:),curveAlg(2,:),'linewidth',2,'Color', colors(2,:)); 
title(sprintf('erf fits\nhuman: %0.02f alg: %0.2f',thresholdHum,thresholdAlg)); 

%% Plot the raw data traces - artifact estimated by the human
centerChannel = elecResp.cells.recElec;
movieIndexRef = find(elecResp.analysis.successRates>0.2,1,'first');
for ii = 1:4
    movieIndex = movieIndexRef - 2 + ii;
    if movieIndex <= length(elecResp.stimInfo.movieNos)
        % Handle the case where stimulation amplitudes are saved as
        % separate movie files.
        if length(unique(elecResp.stimInfo.stimAmps)) == length(elecResp.stimInfo.stimAmps) % There are not 'movies' with duplicate stimulation amps
            subplot(2, 5, 1+ii);
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
            neuronTemplate = elecRespAuto.neuronInfo.templates{1};
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
            subplot(2, 5, 1+ii);
            % Read in all the info from the particular stimulation amp
            if ii == 1
                adjustedIndex = find(elecResp.stimInfo.stimAmps == elecResp.stimInfo.stimAmps(movieIndex));
            else
                adjustedIndex = adjustedIndex + 2;
            end
  
            for aa = 1:length(adjustedIndex)
                dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
                    elecResp.stimInfo.movieNos(adjustedIndex(aa)), 99999);
                subtractionVector = elecResp.analysis.estArtifact{adjustedIndex(aa)}(elecResp.cells.goodElecs == centerChannel, :)';
                latencies = [elecResp.analysis.latencies{adjustedIndex(aa)} elecResp.analysis.otherLatencies{adjustedIndex(aa)}];
                dataToPlot = zeros(size(dataTraces, 1), size(dataTraces, 3));
                nPulses = size(dataTraces, 1);
                for i = 1:nPulses
                    dataToPlot(i, :) = squeeze(dataTraces(i, centerChannel, :)) - subtractionVector;
                end
                
                % Plot neuron template.
                neuronTemplate = elecRespAuto.neuronInfo.templates{1};
                alignmentPoint = find(neuronTemplate == min(neuronTemplate));
                hold on;
                for i = 1:elecResp.stimInfo.nPulses(movieIndex)
                    if ~any(latencies(i,:))
                        misses = plot(dataToPlot(i, :), 'Color',[0.8 0.8 0.8]);
                    else
                        misses = plot(1, 'Color',[0.8 0.8 0.8]);
                    end
                end
                for i = 1:elecResp.stimInfo.nPulses(movieIndex)
                    if any(latencies(i,:))
                        shiftedData = circshift(dataToPlot(i, :)', alignmentPoint - round(latencies(i,:)));
                        hits = plot(shiftedData, 'r');
                    else
                        hits = plot(1,'r'); 
                    end
                end
            end

            plot(neuronTemplate,'k','LineWidth',1);
            ylim([min(neuronTemplate)-50 max(neuronTemplate)+50])
            xlim([0 40]);
            title(sprintf('human\n%0.0f&%0.0f: stimAmp %0.3f uA', adjustedIndex(1),...
                adjustedIndex(2),elecResp.stimInfo.stimAmps(adjustedIndex(1))));
        end
    end
end
try legend([hits misses],'spikes', 'no spikes'); 
catch
end
%% Plot the raw data traces - artifact estimated by the algorithm 
 
if ~isempty(movieIndexRef)
for ii = 1:4    
    movieIndex = movieIndexRef - 2 + ii;
    if size(elecRespAuto.tracesInfo.data,1) == size(elecResp.stimInfo.movieNos,2)/2
        subplot(2, 5, 6+ii);
        adjustedIndex = round(movieIndexRef/2) - 2 + ii; 
        dataTraces = elecRespAuto.tracesInfo.data{adjustedIndex};
        subtractionVector = elecRespAuto.Artifact{1}(adjustedIndex,:);
        dataToPlot = zeros(size(dataTraces, 1), length(subtractionVector));
        nPulses = size(dataTraces, 1);
        for i = 1:nPulses
            dataToPlot(i, :) = squeeze(dataTraces(i, 1:length(subtractionVector))) - subtractionVector;
        end
        latencies = elecRespAuto.latencies{1}(adjustedIndex,:)'; %[elecResp.analysis.latencies{movieIndex} elecResp.analysis.otherLatencies{movieIndex}];
        
        % Plot neuron template.
        neuronTemplate = elecRespAuto.neuronInfo.templates{1};
       
        alignmentPoint = find(neuronTemplate == min(neuronTemplate));
        hold on;
        for i = 1:elecResp.stimInfo.nPulses(adjustedIndex)-1
            if ~any(latencies(i,:))
                misses = plot(dataToPlot(i, :), 'Color',[0.8 0.8 0.8]);
            end
        end
        for i = 1:elecResp.stimInfo.nPulses(adjustedIndex)-1
            if any(latencies(i,:))
                shiftedData = circshift(dataToPlot(i, :)', alignmentPoint + 1 - round(latencies(i,:)));
                hits = plot(shiftedData, 'r');
            end
        end
        plot(neuronTemplate,'k','LineWidth',1);
        ylim([min(neuronTemplate)-50 max(neuronTemplate)+50])
        xlim([0 40]);
        title(sprintf('algorithm\n%0.0f: stimAmp %0.3f uA', adjustedIndex,...
            elecRespAuto.stimInfo.listAmps(adjustedIndex,1)));
    elseif movieIndex <= length(elecRespAuto.tracesInfo.data)
        subplot(2, 5, 6+ii);
        
        dataTraces = elecRespAuto.tracesInfo.data{movieIndex}; 
        if size(elecRespAuto.tracesInfo.data,1) ~= size(elecResp.stimInfo.movieNos,2)
            keyboard; %test these cases, 
        end
        subtractionVector = elecRespAuto.Artifact{1}(movieIndex,:);
        dataToPlot = zeros(size(dataTraces, 1), length(subtractionVector));
        nPulses = size(dataTraces, 1);
        for i = 1:nPulses
            dataToPlot(i, :) = squeeze(dataTraces(i, 1:length(subtractionVector))) - subtractionVector;
        end
        latencies = elecRespAuto.latencies{1}(movieIndex,:)'; %[elecResp.analysis.latencies{movieIndex} elecResp.analysis.otherLatencies{movieIndex}];
        
        % Plot neuron template.
        neuronTemplate = elecRespAuto.neuronInfo.templates{1};
        alignmentPoint = find(neuronTemplate == min(neuronTemplate));
        hold on;
        for i = 1:elecResp.stimInfo.nPulses(movieIndex)-1
            if ~any(latencies(i,:))
                misses = plot(dataToPlot(i, :), 'Color',[0.8 0.8 0.8]);
            end
        end
        for i = 1:elecResp.stimInfo.nPulses(movieIndex)-1
            if any(latencies(i,:))
                shiftedData = circshift(dataToPlot(i, :)', alignmentPoint + 1 - round(latencies(i,:)));
                hits = plot(shiftedData, 'r');
            end
        end
        plot(neuronTemplate,'k','LineWidth',1);
        ylim([min(neuronTemplate)-50 max(neuronTemplate)+50])
        xlim([0 40]);
        title(sprintf('algorithm\n%0.0f: stimAmp %0.3f uA', movieIndex,...
            elecResp.stimInfo.stimAmps(movieIndex)));
    end
end
end
