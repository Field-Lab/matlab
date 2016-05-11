function [eistart eiend] = get_axon_ei_times(ei, threshold)
% usage: [eistart eiend] = get_ei_times(ei)

% Global maximum assumed to be just before start of axon part of EI.
[~, eistart] = max(max(ei));

% EI minimum at each time point
eimin = min(ei);

% From the global minimum, propagate forward as long as eimin is below
% threshold.
for i = eistart:length(eimin)
    if eimin(i) > threshold, break; end
end
eiend = i;
