function [inds indrgbs] = parse_single_cone_rgb_indices(rgbs)
% PARSE_SINGLE_CONE_RGB_INDICES
% usage: [inds indrgbs] = parse_single_cone_rgb_indices(rgbs)
%
% Find stimuli where only a single cone was illuminated
%
% 2012-02 phli
%

inds = {};
for i = 1:numel(rgbs)
    rgb = rgbs{i};
    lum = sum(rgb,2);
    blank = lum == 0;
    unblank = ~blank;
    
    % We only want stimuli where a single cone was illuminated
    if sum(unblank) ~= 1, continue; end
    
    thiscone = find(unblank);
    if length(inds) < thiscone
        inds{thiscone} = [];
        indrgbs{thiscone} = [];
    end
    inds{thiscone}(end+1) = i;
    indrgbs{thiscone}(end+1,:) = rgb(thiscone,:);
end

inds = inds';
indrgbs = indrgbs';