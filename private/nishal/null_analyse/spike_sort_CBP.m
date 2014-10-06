waveforms=cell(9,1);
spike_times=cell(9,1);
spike_amps=cell(9,1);
recon_snippets=cell(9,1);
no_cells=input('Number of cells');
noise=input('Noise Thresold');
for imov=3:9
data=allElecData{imov};
dt=1/20000;
samplingRate=20000;
close all
save('/Volumes/Analysis/nishal/null_data.mat','data','dt');
cd '../CBPSpikesortDemo-master/spikesort_demo/'
[waveforms{imov},spike_times{imov},spike_amps{imov},recon_snippets{imov}]=analNullSpikeSort('null_data',no_cells,noise); % need to give number of potential cells and noise amount .. 
% Store spike sorting result ? 
cd ('../../null_analyse');
save(sprintf('/Volumes/Analysis/nishal/cell%d_2.mat',vision_id),'waveforms','spike_times','spike_amps','recon_snippets','trialStartIndices');
end