

clear dataset

% parameters

% raw data
dataset.raw_spikes = [(1:size(spikeVectors,1))'  spikeVectors];

% num samples
dataset.samples_before = 10;
dataset.samples_after = 14;
dataset.electrodes = 1:7;


% fill in other parameters
dataset.window_length = dataset.samples_before + dataset.samples_after + 1;
dataset.triggers = [];
dataset.mdf_file = [];



[dataset] = spike_plot_and_select_spikes(dataset,'pca',1:7,2,1);

proj_struct = dataset.pca{1};
proj_struct_new = spike_select_several_clusters(proj_struct,5,3);

