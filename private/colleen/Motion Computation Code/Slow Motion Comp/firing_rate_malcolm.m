data_set = '2005-04-26-0';
run = '0';
cell_type = 'off parasol';
cell_indices = get_cell_indices(datarun, {cell_type});
N = length(cell_indices);
T = 900; % length of run
firing_rate = zeros(N,1);

for i=1:N
    firing_rate(i) = length(datarun.spikes{cell_indices(i)})/T;
end

histfit(firing_rate);
title(sprintf('%s run %s %s firing rate',data_set,run,cell_type));
xlabel('Hz');