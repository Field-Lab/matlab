function ipiSpikes = ipiSpikeSort(data, template, recElec, spikeDetThresh, varargin)
%
% quick and dirty spike sorting, used do detect (large-signal) responses in
% intervals between pulses of spatiotemporal pulse sequences (e.g. moving
% bar stimulus)
%


p = inputParser();

p.addRequired('data', @isnumeric) %64 channels x nSamples
p.addRequired('template', @isnumeric) %ei, 64 channels x nSamples
p.addRequired('recElec', @isnumeric) %channel with good signal (on which template was recorded)
p.addRequired('spikeDetThresh', @isnumeric) %thresholds for detecting spikes, as a fraction of the peak ei (template) signal


p.addParamValue('filterTemplate', true, @islogical) %whether or not to high-pass filter the template so that it matches spike waveforms in IPI data that has been HP filtered

p.parse(data, template, recElec, spikeDetThresh, varargin{:})


%% high-pass filter eis to match interpulse spike waveforms 
%(interpulse data filtered using same HP filter)

% if size(template,2) == 1
%     template = template';
% end

if p.Results.filterTemplate;
    c = kaiserord([100 800],[0 1], [0.001 0.001], 20000, 'cell');
    b = fir1(c{:});
    shiftSize = floor(length(b)/2); %how much the trace is shifted by the filtering operation (to correct)
    
    templatePadded = [template zeros(size(template,1), shiftSize*2)];
    
    tmp = filter(b,1,templatePadded,[],2);
    template = tmp(:,shiftSize+1:shiftSize+length(template));
end

%% determine spike times

shiftStep = 0.25;
%for ii = 1:length(recElecs)

spikeDetThreshDAQ = spikeDetThresh*abs(min(template(recElec,:)));

%pad end of data with zeros to avoid index out-of-bounds error
%(this is a hack and hasn't been thoroughly vetted)
data(:,end+1:end+20) = 0;

threshCross = data(recElec,:) < -spikeDetThreshDAQ;
threshCrossStarts = find(diff(threshCross) == 1);

%find location of trace minimum for each threshold crossing event
%(checks over next 20 samples)
threshCrossMins = zeros(size(threshCrossStarts));
for jj = 1:length(threshCrossStarts)
    minVal = min(data(recElec,threshCrossStarts(jj):threshCrossStarts(jj)+20));
    threshCrossMins(jj) = (threshCrossStarts(jj)-1)+find(data(recElec,threshCrossStarts(jj):threshCrossStarts(jj)+20)==minVal,1);
end

% get rid of duplicates
threshCrossMins(diff(threshCrossMins)==0) = []; 

%figure
%hold on

templateMinPos = find(template(recElec,:) == min(template(recElec,:)));
traces = zeros(length(threshCrossMins),31);
tracesAllElec = zeros(64,length(threshCrossMins),31);

%range of offsets = +/- 5 samples from the expected alignment
shiftStart = 11-templateMinPos-5;
shiftEnd = 11-templateMinPos+5;

minErrSpikeTime = zeros(1,length(threshCrossMins));
for jj = 1:length(threshCrossMins)
    fprintf('%c', '.')
    if mod(jj,100)==0
        fprintf('\n')
    end
    traces(jj,:) = data(recElec,threshCrossMins(jj)-10:threshCrossMins(jj)+20);
    zeroRegion = ~any(data(:,threshCrossMins(jj)-10:threshCrossMins(jj)+20),1);
    
    [subtractedTraces offsets] = subtractWithShifts(traces(jj,:), {template(recElec,:)}, shiftStart, shiftEnd, shiftStep);
    subtractedTraces = cell2mat(subtractedTraces);
    subtractedTraces(1,:) = []; offsets(1) = [];
        
    %remove regions of residual corresponding to zero-ed out parts of data
    subtractedTraces(:,zeroRegion) = 0;
    errors = sum(subtractedTraces.^2, 2);
    minErrOffset = offsets(find(errors == min(errors),1));
    minErrSpikeTime(jj) = threshCrossMins(jj)-11+minErrOffset + templateMinPos-1;
    
    %extract new traces based on updated spike times
    %traces(jj,:) = data(recElec,(-10:20) + round(minErrSpikeTime(jj)));
    
    if round(minErrSpikeTime(jj))-10 <= 0
        dataPadded = [zeros(size(data,1),10) data];
        tracesAllElec(:,jj,:) = dataPadded(:,(0:30) + round(minErrSpikeTime(jj)));
        zeroTraceRegions(jj,:) = ~any(dataPadded(:,(0:30) + round(minErrSpikeTime(jj))),1);
    else
        tracesAllElec(:,jj,:) = data(:,(-10:20) + round(minErrSpikeTime(jj)));
        zeroTraceRegions(jj,:) = ~any(data(:,(-10:20) + round(minErrSpikeTime(jj))),1);
    end
end

clear traces

fprintf('\n')

%iteratively pick off the true spikes
realSpikesBin = false(1,size(tracesAllElec,2));
selected = chooseTracesGui(squeeze(tracesAllElec(recElec,:,:)), 'plotAlso', template(recElec,templateMinPos-11:templateMinPos+19),...
    'annotation', ['keep traces that you can tell are real spikes and remove everything else.'...
        10 'If no real spikes are left, remove all traces to exit gui']);
while ~isempty(selected)
    realSpikesBin(selected) = true;
    if all(realSpikesBin) %no traces left
        break
    end
    remainingTraces = reshape(tracesAllElec(recElec,~realSpikesBin,:),sum(~realSpikesBin),[]);
    remainingTraceInds = find(~realSpikesBin);
    selectedTmp = chooseTracesGui(remainingTraces, 'plotAlso', template(recElec,templateMinPos-11:templateMinPos+19),...
        'annotation', ['keep traces that you can tell are real spikes and remove everything else.'...
        10 'If no real spikes are left, remove all traces to exit gui']);
    selected = remainingTraceInds(selectedTmp);
end
clear remainingTraces remainingTraceInds selectedTmp selected

realSpikes = find(realSpikesBin);

%check to make sure none slipped by and no incorrectly identified spikes
%are present
clusterElecs = getCluster(recElec);
figure('position', [100 300 60+length(clusterElecs)*130 300])
yLimAll = [min(template(recElec,templateMinPos-11:templateMinPos+19))...
    max(template(recElec,templateMinPos-11:templateMinPos+19))];

for ii = 1:length(clusterElecs)
    thisElec = clusterElecs(ii);
    h = axes('units', 'pixels', 'position', [130*(ii-1)+30 160 100 100]);
    set(h, 'units', 'normalized')
    
    hold on
    plot(template(thisElec,templateMinPos-11:templateMinPos+19), 'linewidth', 2, 'color', [0.8 0.8 0.8])
    for jj = 1:length(realSpikes)
        x = 1:length(zeroTraceRegions(realSpikes(jj),:));
        plot(x(~zeroTraceRegions(realSpikes(jj),:)), squeeze(tracesAllElec(thisElec, realSpikes(jj), ~zeroTraceRegions(realSpikes(jj),:))),'r')
    end
%     jj = 14
%     if length(realSpikes) >= jj        
%         x = 1:length(zeroTraceRegions(realSpikes(jj),:));
%         plot(x(~zeroTraceRegions(realSpikes(jj),:)), squeeze(tracesAllElec(thisElec, realSpikes(jj), ~zeroTraceRegions(realSpikes(jj),:))),'b')
%         %keyboard
%     end
    
    title('spikes')

    set(h, 'ylim', yLimAll)
    
    axes('units', 'pixels', 'position', [130*(ii-1)+30 30 100 100]); hold on
    set(gca, 'units', 'normalized')
    plot(template(thisElec, templateMinPos-11:templateMinPos+19), 'linewidth', 2, 'color', [0.8 0.8 0.8])
    plot(squeeze(tracesAllElec(thisElec,~realSpikesBin,:))','k')
    title('junk')
    set(gca, 'ylim', yLimAll)
end

%separate plot for incomplete spike waveforms
iInc = 0;
figure('visible', 'off'); hold on
for jj = 1:length(realSpikes)
    x = 1:length(zeroTraceRegions(realSpikes(jj),:));
    if any(zeroTraceRegions(realSpikes(jj),:));
        if zeroTraceRegions(realSpikes(jj),1) %starts with zeroed region
            edgePos = find(~zeroTraceRegions(realSpikes(jj),:),1,'first');
        else %ends with zeroed region
            edgePos = find(~zeroTraceRegions(realSpikes(jj),:),1,'last');
        end
        %if edgePos >= 5 && edgePos <= 15
        if edgePos >= 5 && edgePos <= 15 && any(zeroTraceRegions(realSpikes(jj),11:13))%edge is in region of spike waveform and zeroed region includes spike min position
            iInc = iInc+1;
            plot(100*iInc+template(recElec, templateMinPos-11:templateMinPos+19), 'linewidth', 2, 'color', [0.8 0.8 0.8])
            plot(x(edgePos), 100*iInc+tracesAllElec(recElec,realSpikes(jj),edgePos),'ko')
            plot(x(~zeroTraceRegions(realSpikes(jj),:)),...
                100*iInc+squeeze(tracesAllElec(recElec,realSpikes(jj),~zeroTraceRegions(realSpikes(jj),:))),'r')
        end
    end
end
if iInc > 0
    set(gcf, 'position', [100 100 300 min([iInc*100 800])], 'visible', 'on')
else
    close(gcf)
end

ipiSpikes.spikeTimes = minErrSpikeTime(realSpikesBin);
ipiSpikes.spikeDetThreshDAQ = spikeDetThreshDAQ;




