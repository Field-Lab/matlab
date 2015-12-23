startup_null_analyse_tenessee
%%


datafile = '2014-08-20-0/data001';
type_name= cell(1,1);
type_name{1}='On';

datarun=load_data(datafile)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);
%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!

matlab_cell_ids=get_cell_indices(datarun,type_name);
stas=datarun.stas.stas(matlab_cell_ids);
n_cell=length(stas);


%% Load STAS
% 
% stas_new=cell(length(stas),1);
% for icell=1:length(stas)
%     st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
%     for itime=1:30
%         st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
%     end
%     stas_new{icell}=st_temp;
% end
% stas=stas_new;
% filt_len=size(stas{1},4);

%% 
vision_id=752;
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

init_waveform_all_channels=datarun.ei.eis{matlab_id}(elecList,:);
figure;
plot(init_waveform_all_channels');
%% Load more data

% raster of data
dataTrials=cell(9,1);
allElecData=cell(9,1);
trialStartIndices=cell(9,1);

for imov=[3:9]
    imov
%dataTrials{imov}=getTrialData(sprintf('/Volumes/Data/2014-08-20-0/data00%d',imov),center_elec,330,30);
[dataTrials{imov},allElecData{imov},trialStartIndices{imov}]=getTrialDataMultipleElectrodes(sprintf('/Volumes/Data/2014-08-20-0/data00%d',imov),center_elec,330,30,elecList);
end

  %%
  figure;
 plot(dataTrials{3}(1,800:8000));
 spikeLowerThreshold=input('Spike Lower Threshold');
  
  binSize=40;
 
for imov=[3:9]
   imov
  [~,dataTrialsSpike{imov}]=getTrialSpikes(dataTrials{imov}(:,1:200000),spikeLowerThreshold,'min',binSize);
end

figure; 
for imov=3:9
subplot(7,1,imov-2);
plotSpikeRaster(logical(dataTrialsSpike{imov}));
title(sprintf('Movie %d',imov));
end

%% PSTH
psthBinSize=20;
psthData=cell(9,1);
timeLog=cell(9,1);
for imov=3:9
[timeLog{imov},psthData{imov}]=  psth_calc(dataTrialsSpike{imov},psthBinSize,'nonoverlap');%mean(dataTrialsSpike{imov},1);
end


figure;
for imov=3:9
subplot(7,1,imov-2);
plot(timeLog{imov}*binSize/(20000),psthData{imov});
ylim([0,200]);
title(sprintf('Movie %d .. %0.4f',imov,var(psthData{imov})));
end

pltIdx=[1:length(psthData{3})];
figure;
plot(timeLog{3}*binSize/(20000),psthData{3}(pltIdx),'b')
hold on
plot(timeLog{4}*binSize/(20000),psthData{4}(pltIdx),'r')
hold on
plot(timeLog{6}*binSize/(20000),psthData{6}(pltIdx),'g')
legend('3','4','6')


figure;
plot(timeLog{7}*binSize/(20000),psthData{7}(pltIdx),'b')
hold on
plot(timeLog{8}*binSize/(20000),psthData{8}(pltIdx),'r')
hold on
plot(timeLog{9}*binSize/(20000),psthData{9}(pltIdx),'g')
legend('7','8','9')


%% Do spike sorting? 
spike_sort_CBP
%
%% From STAs, calculate expected firing rate profile.

  %% Make raster
