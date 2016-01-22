function [stimulus, spikes] = block_data_prep(datarun, stimulus_name, block_frames, varargin)
% block frames is a vector of the different block lengths, in order
% = [1200 1200*3] a

p = inputParser;
p.addParameter('block order', 0)
p.addParameter('cell_ids', 0)
p.parse

end