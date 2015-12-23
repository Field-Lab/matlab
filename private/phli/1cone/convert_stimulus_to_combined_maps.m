function newstim = convert_stimulus_to_combined_maps(oldstim)
% CONVERT_STIMULUS_TO_COMBINED_MAPS     Combine maps and update all fields
%
% Old style stimuli were multiple map files with a single index.  New style
% is a single map file with multiple indices.  This updates the maps from
% the old style to the new and then corrects the other stimulus fields to
% match.
%
% The raw pulses are not updated, only the more processed fields, therefore
% raw pulses are dropped from the stuct to avoid confusion.
%
% 2013-03 phli
%

% Simple copy over
newstim.spec = oldstim.spec;

% Combine maps, recalculate nyc and poly versions
newstim.mapims = {combine_maps(oldstim.mapims)};
newstim.mapnyc     = {cellfun(@(C)(C{1}), oldstim.mapnyc,     'UniformOutput', false)'};
newstim.mapnycpoly = {cellfun(@(C)(C{1}), oldstim.mapnycpoly, 'UniformOutput', false)};

% Combine other processed fields
newstim.numcones = length(oldstim.mapindices);
maps = oldstim.maps;
rgbs = oldstim.rgbs;

% Put RGBs in cell array format.  I think the issue is whether it was
% parsed as pulse-sequence or pulse stimulus type...
if ~iscell(rgbs)
    rgbs = num2cell(rgbs, 2)';
end

newstim.rgbs = cell(size(rgbs));
for i = 1:length(newstim.rgbs)
    newstim.rgbs{i} = zeros(newstim.numcones,3);
    newstim.rgbs{i}(maps(i)+1,:) = rgbs{i};
end