function combined = combine_maps(maps)
% COMBINE_MAPS    Combine multiple indexed maps together
% usage: combined = combine_maps(maps)
%
% This handles the special case where there are multiple maps each with
% one index.  We used this in some cases like allcones where we had not
% updated the stimulus map code yet to use a single map with multiple
% indices.
%
% 2012-07 phli
%
combined = maps{1};
for i = 2:length(maps)
    combined = combined + i.*maps{i};
end