clear all

arrayName = 'KM#32';
platType = 'CC';

pathToRawData = '/Data.noindex/Lauren/array_testing/2010-12-22/KM32/data005/data005000.bin';
startSample = 0;
nSamples = 600000; %length of datarun in samples
pulseWidth = 200; %length of test pulse in samples

threshold = 15; %DAQs above mean to trigger at



rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(pathToRawData);

rawData=int16(rawFile.getData(startSample, nSamples))'; %load data: channels x samples

rawData = rawData(2:65,:); %throw out channel 0 data

%% extract mean signal on each electrode

signal = zeros(64, pulseWidth + 120);

for ii = 1:64
    if ~any([9 25 57] == ii)
        aboveThresh = (rawData(ii,:) - mean(rawData(ii,:))) > threshold;
        triggerLoc = find(diff(aboveThresh) == 1) + 1;
        triggersUnused = [];
        

        for jj = length(triggerLoc):-1:1 %remove triggers from noisy downward sloping regions or too close to end of run
            if triggerLoc(jj)+pulseWidth+99 > nSamples ||...
                    sum(aboveThresh(triggerLoc(jj)+1:triggerLoc(jj)+pulseWidth/2)) ~= pulseWidth/2 ||...
                    sum(aboveThresh(triggerLoc(jj)-30:triggerLoc(jj)-1)) ~= 0
%                 if triggerLoc(jj)+pulseWidth+99 <= nSamples
%                     plot(rawData(ii,triggerLoc(jj)-20:triggerLoc(jj)+pulseWidth+99) - mean(rawData(ii,:)));
%                     keyboard
%                 end

                triggersUnused = [triggersUnused triggerLoc(jj)];
                triggerLoc(jj) = [];
                %close

%             else
%                 if triggerLoc(jj)+pulseWidth+99 <= nSamples
%                     plot(rawData(ii,triggerLoc(jj)-20:triggerLoc(jj)+pulseWidth+99) - mean(rawData(ii,:)));
%                     keyboard
%                 end
             end
        end
        
        %hack: throw out first trigger
        triggerLoc = triggerLoc(2:end);
        
        if ii == 1
            triggerLocSave = triggerLoc;
        end
        
        if ii == 10
            triggerLoc = triggerLocSave;
        end
        
        signalSum = zeros(1, pulseWidth+120);
        if 1 %plot to check whether triggers were correctly identified
            figure;
            subplot(1,2,1); hold on
            for jj = 1:length(triggerLoc)
                plot(rawData(ii,triggerLoc(jj)-20:triggerLoc(jj)+pulseWidth+99));
                signalSum = signalSum + double(rawData(ii,triggerLoc(jj)-20:triggerLoc(jj)+pulseWidth+99));
            end
            title(['electrode ' num2str(ii)])
            hold off

            subplot(1,2,2); hold on
            for jj = 1:length(triggersUnused)
                if triggersUnused(jj)+pulseWidth+99 > nSamples
                    plot(rawData(ii,triggersUnused(jj)-20:end),'k');
                elseif triggersUnused(jj)-20 < 1
                    plot(-triggersUnused(jj)+22:(pulseWidth+120), rawData(1:(triggersUnused(jj)+pulseWidth+99)),'k')
                else
                    plot(rawData(ii,triggersUnused(jj)-20:triggersUnused(jj)+pulseWidth+99),'k');
                end
            end

            hold off
        end
        
        
        signal(ii,:) = signalSum/length(triggerLoc);
    end
end

%% generate plots

[xCoords yCoords] = getElectrodeCoords61();

arrayWidth = max(xCoords)*2*1.2;
arrayHeight = max(yCoords)*2*1.2;

figure('position', [100 100 600 525])

yMax = max(max(signal(signal ~= 0)));
yMin = min(min(signal(signal ~= 0)));

maxRange = 0;
ranges = zeros(1,64);
for ii = 1:64
    ranges(ii) = max(signal(ii,:))-min(signal(ii,:));
    maxRange = max(maxRange, ranges(ii));
end

rangesNorm = ranges/maxRange;

for i = 1:64
    if ~(i==9||i==25||i==57)
        axes('position', [xCoords(i)/arrayWidth + 0.45, yCoords(i)/arrayHeight + 0.45, 0.1, 0.1])

        plot(signal(i,:), 'k')
        set(gca, 'ylim', [yMin yMax])
        %set(gca, 'ylim', [-10 8], 'xtick', [], 'ytick', [])
        axis off
    end
end

axes('position', [0.05 0.9 0.3 0.1])
text(0, 0.5, [arrayName, 10, 'platinized with ' platType, 10 'max signal range = ' num2str(maxRange)])
axis off

% array orientation legend
axes('position', [0.8 0.8 0.15 0.15])
hold on
plot(xCoords([8 19 30 40 51 62 8]), yCoords([8 19 30 40 51 62 8]), 'k-')
text(xCoords(7),  yCoords(7),  '8',  'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(61), yCoords(61), '62', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(50), yCoords(50), '51', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(39), yCoords(39), '40', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(29), yCoords(29), '30', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(18), yCoords(18), '19', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
hold off
axis equal
axis off


%% figure using dots only

figure('position', [100 100 600 525])
hold on
axis equal

for i = 1:64
    if ~(i==9||i==25||i==57)
        if rangesNorm(i)>.02
            plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(rangesNorm(i)*40), 'MarkerFaceColor', 'k')
        end
    end
end
axis off

axes('position', [0.05 0.9 0.3 0.1])
text(0, 0.5, [arrayName, 10, 'platinized with ' platType, 10 'max signal range = ' num2str(maxRange)])
axis off


% array orientation legend
axes('position', [0.8 0.8 0.15 0.15])
hold on
plot(xCoords([8 19 30 40 51 62 8]), yCoords([8 19 30 40 51 62 8]), 'k-')
text(xCoords(7),  yCoords(7),  '8',  'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(61), yCoords(61), '62', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(50), yCoords(50), '51', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(39), yCoords(39), '40', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(29), yCoords(29), '30', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(18), yCoords(18), '19', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
hold off
axis equal
axis off




