% rawDataDir = '/Volumes/Acquisition/Data/2014-12-10-SB1test/';
% datarun = 'data000'; 
rawDataDir = '/Volumes/Acquisition/Data/2014-12-10-noisetests/';
datarun = 'data010'; 
rawDataPath = [rawDataDir datarun];
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(rawDataPath);
header = rawFile.getHeader();
totalSamples = header.getNumberOfSamples(); 
% channels = 1:(header.getNumberOfElectrodes-1); %-1 because first channel has TTL pulses
% NS_GlobalConstants = NS_GenerateGlobalConstants(length(channels));

%%
plottingOn = 1;
Fs = 20000;                   % Sampling frequency
T = 1/Fs;                     % Sample time
L = 10000;                    % Length of signal
t = (0:L-1)*T;                % Time vector


%% Reading in the data 
startSample = 0;
nSamplesToRead = totalSamples-startSample;
%%
figure;

while nSamplesToRead>0
    tic
    
    rawData = double(rawFile.getData(startSample,min(10000,nSamplesToRead)));
    nSamplesToRead = nSamplesToRead - 10000;
    startSample = startSample + 10000;
    perChannelThresholds = median(rawData) - 30; 
    [rowIdx,colIdx] = find(rawData < repmat(perChannelThresholds,size(rawData,1),1));
    cla; 
    plot(1000*t(1:2:10000),rawData(1:2:10000,2:end)); 
    pause(0.1); 
%     if colIdx
%         plot(1000*t(1:2:10000),rawData(1:2:10000,colIdx))
%         %     plot(1000*t(1:2:10000),rawData(1:2:10000,2:end))
%         title(sprintf('raw signal from all below threshold, %d channels',...
%             size(unique(colIdx),1)));
%         xlabel('time (milliseconds)');
%         pause(0.01);
%     end
    fprintf('\n%d samples to read; %0.1f%% complete\n',nSamplesToRead,100*startSample/totalSamples)
    toc
%     rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(rawDataPath);
%     header = rawFile.getHeader();
%     totalSamples = header.getNumberOfSamples();
if size(unique(colIdx),1)>10
    disp('Bad run');
end
end
rawFile.close; 
%% 
positions = loadElecPositions512()
perChannelThresholds = median(rawData) - 30; 
[I,J] = find(rawData < repmat(perChannelThresholds,size(rawData,1),1))
figure; plot(1000*t(1:2:10000),rawData(1:2:10000,J))
    title('raw signal from all channels')
    xlabel('time (milliseconds)');
    title(sprintf('%d samples to read',nSamplesToRead));
figure; scatter(positions(:,1),positions(:,2),100,...
    abs(perChannelThresholds(:,2:end)),'filled'); 
%% testing stuff
figure;
for chan = 1:512
    cla;
    plot(rawData(:,1+chan)); hold on; 
    line([0 10000],[median(rawData(:,1+chan)) median(rawData(:,1+chan))],...
        'Color','red'); title(num2str(chan))
    line([0 10000],[median(rawData(:,1+chan))-30 median(rawData(:,1+chan))-30],...
        'Color','green'); title(num2str(chan));
    ylim([-400 -200]);
    pause(0.1); 
end