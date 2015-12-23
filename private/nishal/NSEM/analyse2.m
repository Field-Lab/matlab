startup_analyse_tenessee
%%


datafile = '2012-08-09-3/data005-from-data002/data005-from-data002';
imov=5;
bin_datafile=sprintf('/Volumes/Data/2012-08-09-3/data00%d',imov);
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



    imov
%dataTrials{imov}=getTrialData(sprintf('/Volumes/Data/2014-08-20-0/data00%d',imov),center_elec,330,30);
[allElecData{imov}]=getTrialDataMultipleElectrodesOnly(bin_datafile,center_elec,1.5*60*60,elecList);





%% Do CBP spike sorting? 
waveforms=cell(9,1);
spike_times=cell(9,1);
spike_amps=cell(9,1);
recon_snippets=cell(9,1);
no_cells=input('Number of cells');
noise=input('Noise Thresold');
data=allElecData{imov};
dt=1/20000;
samplingRate=20000;
close all
save('/Volumes/Analysis/nishal/nsem_data.mat','data','dt','-v7.3');
cd '/Volumes/Analysis/nishal/CBPSpikesortDemo-master/spikesort_demo/'
[waveforms{imov},spike_times{imov},spike_amps{imov},recon_snippets{imov}]=nsemSpikeSort('nsem_data',no_cells,noise); % need to give number of potential cells and noise amount .. 
% Store spike sorting result ? 
cd ('../../NSEM');
save(sprintf('/Volumes/Analysis/nishal/NSEM_cell%d_long.mat',vision_id),'waveforms','spike_times','spike_amps','recon_snippets');
save(sprintf('/Volumes/Analysis/nishal/NSEM_cell%d_long2.mat',vision_id),'-v7.3');

%% 
figure;
for icell=1:no_cells
subplot(3,2,icell);
plot(waveforms{imov}{icell});
end

subplot(3,2,6);
plot(init_waveform_all_channels');

correct_cell=input('Correct cell ? ');
spk_tm=spike_times{imov}{correct_cell};
spk_amp=spike_amps{imov}{correct_cell};
save(sprintf('/Volumes/Analysis/NSEM_cell_%d_long_share.mat',vision_id),'spk_tm','spk_amp');


