
cd /Volumes/Analysis/nishal/NSEM_SS/

load('../CBPcells.mat');
load('NSEM_dataset_2013-10-10-0_cell346_data001_a.mat');

no_cells=length(spike_times{imov});

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
save(sprintf('/Volumes/Lab/Users/bhaishahster/SS/SS_dataset_%s_cell%d_data00%d_share.mat',dataset,vision_id,imov),'spk_tm','spk_amp');
