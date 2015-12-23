function spike_frame_offset_hist(datarun, cell_spec, offset_depth, bins, color)
% SPIKE_FRAME_OFFSET_HIST    Make histogram of offset times from frame change to spikes
%
% usage: spike_frame_offset_hist(datarun, cell_spec, offset_depth, bins, color)
%
% 2010-02 phli
%

if nargin < 5
    color = 'b';
end

if nargin < 4
    bins = 500;
end



cell_indices = get_cell_indices(datarun, cell_spec);

if ~isfield(datarun, 'stimulus') || ~isfield(datarun.stimulus, 'frame_times')
    datarun = get_frame_times(datarun);
end
frame_times = datarun.stimulus.frame_times;


for i = 1:numel(cell_indices)
    offsets{i} = calc_spike_frame_offsets(datarun.spikes{cell_indices(i)}, frame_times, offset_depth);
end
offsets = vertcat(offsets{:});


[n,x] = hist(offsets(:), bins);
n = n ./ max(n(:));
%h = bar(x, n, 'hist');
plot(x, n, color);