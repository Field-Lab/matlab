stats_file = '/Volumes/Analysis/2012-08-09-3/data002/data002-stats.txt';
stats = load(stats_file);


color = stats(:,3);
cell_id = stats(:,4);
electrode_number = stats(:,5);
cluster_number = stats(:,6);
duplicate_index = stats(:,7);
firing_rate = stats(:,8);
contamination_index = stats(:,9);
fragmentation_index = stats(:,10);
fragmented_clusters = stats(:,11);
nearest_neighbor = stats(:,12);
best_union = stats(:,13);
worst_likelihood = stats(:,14);
best_likelihood_neighbor = stats(:,15);

