%icaTesterScript

details = getDatasetDetails(1);

nMovies = length(details.movieNumbers);
nPatterns = length(details.patternNumbers);


startFrame = 1;
endFrame = 60;

dataTraces = NS_ReadPreprocessedData(details.dataPath, '', 0, details.patternNumbers(1), details.movieNumbers(10));

nTraces  = size(dataTraces, 1);
nSamples = size(dataTraces, 3);

figure
plot(squeeze(dataTraces(:,details.centerChannel,startFrame:endFrame))')

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);

channelsToUse=electrodeMap.getAdjacentsTo(details.centerChannel, details.channelRadius)';

%% attempt one: repetitions x samples

data = squeeze(dataTraces(1:50, details.centerChannel, startFrame:endFrame));

%data = channels x (frames*epochs)

[weights, sphere] = runica(data);

components = weights*data;

figure
plot(components')
title('using multiple trials')

%plot(data')


%% attempt 2: channels x samples

data = squeeze(dataTraces(1, channelsToUse, startFrame:endFrame));

[weights, sphere] = runica(data);

components = weights*data;

figure
plot(components')
title('using multiple electrodes')

%% attempt 3: channels x multiple epochs

frameLength = endFrame - startFrame + 1;
data = zeros(length(channelsToUse), (frameLength*size(dataTraces,1)));

for i = 1:size(dataTraces,1)
    data(:,(i-1)*frameLength+1 : i*frameLength) = squeeze(dataTraces(i, channelsToUse, startFrame:endFrame));
end

%baseline-zeroing each epoch (remove mean of basevector from each channel and epoch)
[data datamean] = rmbase(data, frameLength, 1:frameLength);

figure
hold on
for i = 1:length(channelsToUse)
    plot(data(i, 1:5*frameLength)+(i-1)*100)
end
hold off
title('first 5 epochs')



[weights, sphere] = runica(data);

mixMat = inv(weights);
components = weights*data;

figure
for i = 1:length(channelsToUse)
    subplot(length(channelsToUse), 1, i)
    plot(components(i,1:10*frameLength))
    if i==1
        title('using multiple electrodes with concatenated epochs')
    end
end


%taking mean over epochs
componentsSplit = shiftdim(reshape(components, length(channelsToUse), frameLength, []),2);
meanComponent = squeeze(mean(componentsSplit));

figure
hold on
for i = 1:length(channelsToUse)
    plot(meanComponent(i, :) + (i-1)*100)
end
hold off
title('mean components')

figure
for i = 1:length(channelsToUse)
    subplot(length(channelsToUse), 1, i)
    plot(meanComponent(i,:))
    if i==1
        title('mean of components using multiple electrodes with concatenated epochs')
    end
end

%% zeroing out some of the components

% mixMat(:,3) = 0;
% mixMat(:,4) = 0;
% mixMat(:,5) = 0;
% mixMat(:,6) = 0;
% mixMat(:,7) = 0;



%% remixing
remixed = mixMat*components;

figure
hold on
for i = 1:length(channelsToUse)
    plot(remixed(i, 1:5*frameLength)+(i-1)*100)
end
hold off
title('remixed')
