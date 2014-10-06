clear all
close all


%% need to set each time

electrodes = [4 6 1 15 41 58];  %the electrodes to be included for PCA

frequency = 10;                 %frequency of administered pulses, in Hz


%% specific to full version (mac)

pathToData = '/Volumes/Lee/Data/Lauren/spike_size_analysis/2008-03-25-5/data001/data001000.bin';

trialLength = 5;               %length of trial, in seconds

%opens data file
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(pathToData);

%gets data from sample: first index is data point, second index is
%electrode number (actual number -1 because of indexing)
data = rawFile.getData(0, 20000*trialLength);
%data = rawFile.getData(0.);

header = rawFile.getHeader();
comment = header.getComment();


%% specific to partial version (pc)

%load analysis\data001.mat data trialLength


%% need to set sometimes

%beginning and end of time window relative to stimulus, in ms
windowStart = 0;
windowEnd = 3;

offset = 25;                    %offset from zero when first pulse occurs, in ms
baselineoffset = 20;            %offset from previous pulse for beginning of baseline voltage measurement for next pulse

disclude = 10;                  %number of initial pulses to be discluded from analysis
            

%% keep the same

melectrodes = electrodes + 1;   %shifted indexing for matlab

totalPulses = frequency*trialLength - disclude;    %number of stimulation pulses during trial

color = cool(totalPulses);

winStartReal = (offset+windowStart)*20 + 1;
winEndReal = (offset+windowEnd)*20;
vectorLength = (windowEnd - windowStart)*20;

samplestemp = zeros(totalPulses, vectorLength);
samplescat = zeros(totalPulses, vectorLength*length(melectrodes));


% places data within specified window of time for each pulse after the
% first 10 into the array "samples"
% *average value of trace between 20 ms after the previous pulse and the
% start of the current pulse is subtracted from data to set baseline to 0

for j = 1:length(melectrodes)
    figure
    hold on
    for i=1:totalPulses
        baselineSample = data(offset*20+baselineoffset*20+1+(i+disclude-2)*20000/frequency : offset*20+1+(i+disclude-1)*20000/frequency, melectrodes(j));
        baseline = mean(baselineSample);
        samplestemp(i,:) = data(winStartReal+(i+disclude-1)*20000/frequency : winEndReal+(i+disclude-1)*20000/frequency, melectrodes(j)) - baseline;
        samplescat(:,vectorLength*(j-1)+1:vectorLength*j) = samplestemp;
        
        current = plot(winStartReal : 1 : winEndReal, samplestemp(i,:));
        set(findobj(current,'Type','line'),'Color',color(i,:))
    end
    hold off
end

[PCAcoef, PCAscore] = princomp(samplescat);

[selx,sely,indexnr] = lasso(PCAscore(:,1), PCAscore(:,2));
%igorindexes = indexnr + 9;

%% PC scatter plots
% 
% figure(
% 
% subplot(2,3,1)
% plot(PCAscore(:,1),PCAscore(:,2),'.')
% xlabel('PC1','FontWeight','bold')
% ylabel('PC2','FontWeight','bold')
% 
% subplot(2,3,2)
% plot(PCAscore(:,1),PCAscore(:,3),'.')
% xlabel('PC1','FontWeight','bold')
% ylabel('PC3','FontWeight','bold')
% 
% subplot(2,3,3)
% plot(PCAscore(:,1),PCAscore(:,4),'.')
% xlabel('PC1','FontWeight','bold')
% ylabel('PC4','FontWeight','bold')
% 
% subplot(2,3,4)
% plot(PCAscore(:,2),PCAscore(:,3),'.')
% xlabel('PC2','FontWeight','bold')
% ylabel('PC3','FontWeight','bold')
% 
% subplot(2,3,5)
% plot(PCAscore(:,2),PCAscore(:,4),'.')
% xlabel('PC2','FontWeight','bold')
% ylabel('PC4','FontWeight','bold')
% 
% subplot(2,3,6)
% plot(PCAscore(:,3),PCAscore(:,4),'.')
% xlabel('PC3','FontWeight','bold')
% ylabel('PC4','FontWeight','bold')
% 
% %% PC weights
% 
% figure(j+2)
% hold on
% plot(PCAcoef(:,1)*2000,'c-')
% plot(PCAcoef(:,2)*2000,'m-')
% plot(PCAcoef(:,2)*2000,'g-')
% plot(PCAcoef(:,4)*2000,'b-')
% for i=1:totalPulses-10
%     plot(samplescat(i,:),'k-')
% end
% plot(PCAcoef(:,1)*2000,'c-','LineWidth',2)
% plot(PCAcoef(:,2)*2000,'m-','LineWidth',2)
% plot(PCAcoef(:,3)*2000,'g-','LineWidth',2)
% plot(PCAcoef(:,4)*2000,'b-','LineWidth',2)
% legend('PC1','PC2','PC3','PC4')
% hold off