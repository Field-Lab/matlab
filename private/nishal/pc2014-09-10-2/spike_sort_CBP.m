waveforms=cell(1,1);
spike_times=cell(1,1);
spike_amps=cell(1,1);
recon_snippets=cell(1,1);
no_cells=input('Number of cells');
noise=input('Noise Thresold');

imov=1;
data=allElecData;
dt=1/20000;
samplingRate=20000;
close all
save('/Volumes/Analysis/nishal/nsem_data.mat','data','dt');
cd '../CBPSpikesortDemo-master/spikesort_demo/'
%[waveforms{imov},spike_times{imov},spike_amps{imov},recon_snippets{imov}]=analNullSpikeSort('nsem_data',no_cells,noise); % need to give number of potential cells and noise amount .. 
[waveforms{imov},spike_times{imov},spike_amps{imov},recon_snippets{imov}]=nsemSpikeSort('nsem_data',no_cells,noise);
% Store spike sorting result ? 
cd ('../../pc2014-09-10-2/');
save(sprintf('/Volumes/Analysis/nishal/cell%d_3.mat',vision_id),'waveforms','spike_times','spike_amps','recon_snippets','init_waveform_all_channels');
%save(sprintf('/Volumes/Analysis/nishal/cell%d.mat',vision_id),'-v7.3');
