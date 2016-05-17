function fitspikes = concat_spikes(spikes_cell, block_length)
fitspikes = [];
for block = 1:length(spikes_cell)
    fitspikes = [fitspikes; spikes_cell{block}+block_length*(block-1)];
end
end