data_set = '2005-04-26-0';
cell_type = 'on on midget nc4';
cell_indices = get_cell_indices(datarun, {cell_type});
tmax = ceil(max(datarun.spikes{1}));

t = 0:0.001:tmax; % 1ms time bins
N = length(cell_indices);
T = 100; % num displacements
spike_autocorr = zeros(N, T+1); 

% iterate over cells
for i = 1:N;
    spikes = hist(datarun.spikes{cell_indices(i)}, t);
    spike_autocorr(i,:) = autocorr(spikes, T);
end

mean_spike_autocorr = mean(spike_autocorr);
plot(1:T,mean_spike_autocorr(2:(T+1)));
title(sprintf('%s %s mean autocorr',data_set,cell_type));
xlabel('ms');