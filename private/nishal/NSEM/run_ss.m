
%%


% type_name= cell(1,1);
% type_name{1}='On Parasol';

datarun=load_data(analysis_datafile)
%datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);

%% 

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



    imov
%dataTrials{imov}=getTrialData(sprintf('/Volumes/Data/2014-08-20-0/data00%d',imov),center_elec,330,30);
[allElecData{imov}]=getTrialDataMultipleElectrodesOnly(bin_datafile,center_elec,elecList);





%% Do CBP spike sorting? 
waveforms=cell(9,1);
spike_times=cell(9,1);
spike_amps=cell(9,1);
recon_snippets=cell(9,1);

data=allElecData{imov};
dt=1/20000;
samplingRate=20000;
close all
save('/Volumes/Analysis/nishal/SS/ss_data.mat','data','dt','-v7.3');
cd '/Volumes/Analysis/nishal/CBPSpikesortDemo-master/spikesort_demo/'
[waveforms{imov},spike_times{imov},spike_amps{imov},recon_snippets{imov}]=nsemSpikeSort('nsem_data',no_cells,noise); % need to give number of potential cells and noise amount .. 
% Store spike sorting result ? 
cd ('~/Nishal/matlab/private/nishal/NSEM');
save(sprintf('/Volumes/Analysis/nishal/SS/SS_dataset_%s_cell%d_data00%d_a_%dcells.mat',dataset,vision_id,imov,no_cells),'waveforms','spike_times','spike_amps','recon_snippets','init_waveform_all_channels','imov','dataset','analysis_datafile','bin_datafile','vision_id','no_cells');
%save(sprintf('/Volumes/Analysis/nishal/NSEM_SS/NSEM_dataset_%s_cell%d_data00%d_b_%dcells.mat',dataset,vision_id,imov,no_cells),'-v7.3');

