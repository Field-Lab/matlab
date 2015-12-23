function [offsets preceeding_frame_indices] = calc_spike_frame_offsets(spikes, frame_times, depth)
% CALC_SPIKE_FRAME_OFFSETS    Determine for each spike how long it followed a frame change by
%
% usage: [offsets preceeding_frame_indices] = calc_spike_frame_offsets(spikes, frame_times, depth)
%
% inputs: SPIKES         The spike times, units compatible with FRAME_TIMES
%         FRAME_TIMES    The frame times, units compatible with SPIKES
%
%         DEPTH          How many frames back to get the offsets from.  By
%                        default, just get offset to immediately preceeding
%                        frame.  However, in general our frame rates are
%                        high enough that the cell is actually responding
%                        to a frame several frames before that.
%
%
% outputs: OFFSETS                     A little tricky.  This will be a NxD
%                                      matrix of the spike offset times.  N
%                                      is the number of spikes, i.e. the
%                                      length of the SPIKES input vector.
%                                      D is the DEPTH, so for each spike
%                                      there will be several offsets
%                                      depending on how many frame changes
%                                      back we were interested in.
%
%          PRECEEDING_FRAME_INDICES    Also NxD.  The indices into
%                                      FRAME_TIMES for the frames that were
%                                      found to preceed each spike.
%
%
% 2010-02 phli
%

if nargin < 3
    depth = 1;
end


% Determine the frames that precede each spike up to DEPTH frames
max_frame = numel(frame_times);
preceeding_frame_indices = zeros(numel(spikes), depth);
frame_index = 1;
for i = 1:numel(spikes)
    while frame_index < max_frame && spikes(i) > frame_times(frame_index)
        frame_index = frame_index + 1;
    end

    preceeding_frame_indices(i,:) = ((frame_index-depth):(frame_index-1))';
end



% Pad frame_times with NaN up to DEPTH; hack to handle the negative indices
frame_times = [repmat(NaN, depth, 1); frame_times(:)];

for i = 1:depth
    offsets(:,i) = spikes - frame_times(preceeding_frame_indices(:,i) + depth); %#ok<AGROW> % AGROW shouldn't matter as depth is probably not more than 15
end
