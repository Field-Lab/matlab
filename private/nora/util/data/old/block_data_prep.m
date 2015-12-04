function  block_data_prep(datarun, block_frames, varargin)
% block frames is a vector of the different block lengths, in order
% = [1200 1200*3] a
% Assumes the first block is a repeated segment and the second block is a fitting
% segment

p = inputParser;
p.addParameter('block_order', 0)
p.addParameter('cell_ids', 0)
p.addParameter('stimulus_name', 0)
p.parse

% define the number of blocks and their order. 
% Usually n_block_types will be 2: fitting and testing
n_block_types = length(block_frames);
n_blocks = 100; % how should this be figured out?
if p.Results.block_order
   block_order = p.Results.block_order;
else
   block_order =  repmat(1:n_block_types,1, ceil(n_blocks/n_block_types));
end

% loading up the appropriate stimulus
if p.Results.stimulus_name
    if strcmp(p.Results.stimulus_name(1:2), 'BW') || strcmp(p.Results.stimulus_name(1:2), 'RG')
        load_movie
    elseif strcmp(p.Results.stimulus_name(1:2), 'NS')
    end
end

% load up the triggers 
triggers = datarun.triggers;
triggers_per_block = block_frames/100;
start_trigger(1) = 1;
for i = 1:n_blocks
    end_trigger = start_trigger(i) + triggers_per_block(block_order(i));
    start_trigger(i+1) = end_trigger + 1;
end

% visual check
plot(diff(triggers))
hold on
plot(start_trigger, 1/1.2*ones(length(start_trigger)), '*')

% start times
start_time = triggers(start_trigger);


end