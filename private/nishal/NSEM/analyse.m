startup_analyse_tenessee
%%


datafile = '2012-08-09-3/data005-from-data002/data005-from-data002';
% type_name= cell(1,1);
% type_name{1}='On Parasol';

datarun=load_data(datafile)
%datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);

%% 
vision_id=2101;
idx=[1:length(datarun.cell_ids)];
matlab_id=idx(datarun.cell_ids==vision_id);
cell_ei=datarun.ei.eis{matlab_id};

cell_ei_mag=sum(cell_ei.^2,2);
[v,center_elec]=max(cell_ei_mag);

% Get all cells which are strong at calculated center elctrode
% otherCellsMatlabID;
magCenterElec=[];
for icell=1:length(datarun.ei.eis)
    x=datarun.ei.eis{icell}(center_elec,:);
    magCenterElec=[magCenterElec;sum(x.^2)];
end
figure;stem(magCenterElec);
xlabel('Matlab Cell ID');
ylabel('Magnitude');
title(sprintf('Magnitude of different cells on center electrode %d',center_elec));

% get surround electrodes 
centerElecPos=datarun.ei.position(center_elec,:);
distCenterElec = datarun.ei.position-repmat(centerElecPos,[size(datarun.ei.position,1),1]);
distCenterElec=sum(distCenterElec.^2,2);
[sortedDist,nearestChannels]=sort(distCenterElec,'ascend');
elecList=nearestChannels(1:7);


% Get Initial waveform estimates
% For center electrode
init_waveform=datarun.ei.eis{matlab_id}(center_elec,:);
figure;plot(init_waveform);
title('Waveform shape on center electrode');

init_waveform_all_channels=datarun.ei.eis{matlab_id}(elecList,:);
figure;
plot(init_waveform_all_channels');
title('Waveform shape on center and 6 surrounding electrodes');

%% Load more data

% raster of data
dataTrials=cell(9,1);
allElecData=cell(9,1);
trialStartIndices=cell(9,1);

for imov=[3]
    imov
%dataTrials{imov}=getTrialData(sprintf('/Volumes/Data/2014-08-20-0/data00%d',imov),center_elec,330,30);
[dataTrials{imov},allElecData{imov},trialStartIndices{imov}]=getTrialDataMultipleElectrodes(sprintf('/Volumes/Data/2012-08-09-3/data00%d',imov),center_elec,1.5*60*60,1,elecList);
end

  %%
  figure;
 plot(dataTrials{3}(1,800:8000));
 title(sprintf('VID %d',vision_id));
 xlabel('Samples')
 ylabel('Voltage on central channel')
 spikeLowerThreshold=input('Spike Lower Threshold');
  
  binSize=40;
 
for imov=[3]
   imov
  [~,dataTrialsSpike{imov}]=getTrialSpikes(dataTrials{imov}(:,1:200000),spikeLowerThreshold,'min',binSize);
end

figure; 
for imov=3:3
subplot(1,1,imov-2);
plotSpikeRaster(logical(dataTrialsSpike{imov}));
title(sprintf('Movie %d',imov));
end

%% PSTH
% psthBinSize=20;
% psthData=cell(9,1);
% timeLog=cell(9,1);
% for imov=3:3
% [timeLog{imov},psthData{imov}]=  psth_calc(dataTrialsSpike{imov},psthBinSize,'nonoverlap');%mean(dataTrialsSpike{imov},1);
% end
% 
% 
% figure;
% for imov=3:3
% subplot(1,1,imov-2);
% plot(timeLog{imov}*binSize/(20000),psthData{imov});
% ylim([0,200]);
% title(sprintf('Movie %d .. %0.4f',imov,var(psthData{imov})));
% end
% 
% pltIdx=[1:length(psthData{3})];
% figure;
% plot(timeLog{3}*binSize/(20000),psthData{3}(pltIdx),'b')
% 

%% Do CBP spike sorting? 
spike_sort_CBP
%
%% Load previous result? 
load('/Volumes/Analysis/nishal/NSEM_cell3676_2.mat');

load('/Volumes/Analysis/nishal/NSEM_cell3676_1.mat');
% have spike_times, spike_amps, etc
%% Take spiking information found from vision

neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile('/Volumes/Analysis/2012-08-09-3/data003-from-data002/data003-from-data002.neurons');
vision_spk_times= neuronFile.getSpikeTimes(vision_id);
  trialStartIndices=trialStartIndices{3};
  %% Make raster

  figure;
  imov=3;
  n_cell=no_cells
  for icell=1:n_cell
      subplot(3,2,icell);
      plot(waveforms{3}{icell});
      title(sprintf('Cell : %d',icell));
  end
  
  subplot(3,2,5);
  plot(init_waveform_all_channels');
  title('Vision Waveform');
  
    correct_cell=input('What is the correct cell?');
   cbp_spk_times=spike_times{imov}{correct_cell};
    
   vision_spk_train =zeros(20000*1100,1);
   vision_spk_train(vision_spk_times)=1;
   
   cbp_spk_train=zeros(20000*1100,1);
   cbp_spk_train(round(cbp_spk_times))=1;
   
      
    dummyTrain=cbp_spk_train;
 %  trialStartIndices=trialStartIndices{3};
    trialIntervalSamples = 11*20000;
    numTrials=100;
    
    % NOTE: use only 29 samples
    dataTrials=zeros(numTrials,trialIntervalSamples);
    for iTrial=1:numTrials-1 % 30th sample less data?
    dataTrials(iTrial,:)=dummyTrain(trialStartIndices(iTrial):trialStartIndices(iTrial)+trialIntervalSamples-1)';
    end

     dataTrials_cbp = dataTrials;
   
       
    dummyTrain=vision_spk_train;
    
    % NOTE: use only 29 samples
    dataTrials=zeros(numTrials,trialIntervalSamples);
    for iTrial=1:numTrials-1 % 30th sample less data?
    dataTrials(iTrial,:)=dummyTrain(trialStartIndices(iTrial):trialStartIndices(iTrial)+trialIntervalSamples-1)';
    end

     dataTrials_vision = dataTrials;
 %%  
binSize=20;

spikesBinnedcbp=getTrialSpikes( dataTrials_cbp,0.01,'max',binSize);
spikesBinnedvision=getTrialSpikes( dataTrials_vision,0.01,'max',binSize);

              LineFormat = struct();
              LineFormat.Color = [0.3 0.3 0.3];
              LineFormat.LineWidth = 0.35;
              LineFormat.LineStyle = ':';
   %           plotSpikeRaster(spikes,'LineFormat',LineFormat)
figure;    
subplot(3,1,1);
plotSpikeRaster(logical(spikesBinnedvision),'LineFormat',LineFormat);
title('Vision');
xlabel('Time');
ylabel('Trial');

subplot(3,1,2);
LineFormat.Color = [1 0 0];
plotSpikeRaster(logical(spikesBinnedcbp),'LineFormat',LineFormat);
title('CBP');
xlabel('Time');
ylabel('Trial');
%% Better Raster

numTrials=100;
dummyTimes=vision_spk_times;
dummyTrialTimes=cell(numTrials-1,1);
for iTrial=1:numTrials-1
dummyTrialTimes{iTrial}=dummyTimes(dummyTimes>=trialStartIndices(iTrial)& dummyTimes<trialStartIndices(iTrial+1))-trialStartIndices(iTrial);
dummyTrialTimes{iTrial}=double(dummyTrialTimes{iTrial}')/20000;
end
vision_spks=dummyTrialTimes;


numTrials=100;
dummyTimes=cbp_spk_times;
dummyTrialTimes=cell(numTrials-1,1);
for iTrial=1:numTrials-1
dummyTrialTimes{iTrial}=dummyTimes(dummyTimes>=trialStartIndices(iTrial)& dummyTimes<trialStartIndices(iTrial+1))-trialStartIndices(iTrial);
dummyTrialTimes{iTrial}=double(dummyTrialTimes{iTrial}')/20000;
end
cbp_spks=dummyTrialTimes;

figure;    
subplot(2,1,1);
plotSpikeRaster(vision_spks);

title(sprintf('Vision , Correct cell: %d, VID: %d',correct_cell,vision_id));
xlabel('Time');
ylabel('Trial');

subplot(2,1,2);
plotSpikeRaster(cbp_spks);
title('CBP');
xlabel('Time');
ylabel('Trial');

%% PSTH
psthBinSize=40%20;

[timeLog,psthData]=  psth_calc(spikesBinnedcbp,psthBinSize,'nonoverlap');%mean(dataTrialsSpike{imov},1);

figure;

plot(timeLog/(20000),psthData);

[timeLog,psthData]=  psth_calc(spikesBinnedvision,psthBinSize,'nonoverlap');%mean(dataTrialsSpike{imov},1);
hold on;
plot(timeLog/(20000),psthData,'r');
legend('CBP','Vision')
title(sprintf('correct cell %d, vision id %d',correct_cell,vision_id));


