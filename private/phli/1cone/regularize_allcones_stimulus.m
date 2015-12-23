function stimulus = regularize_allcones_stimulus(stimulus)
% Get STIMULUS into canonical form, making various assumptions

% Need to parse pulse information into separated maps and rgbs?
if ~isfield(stimulus, 'mapindices')
    stimulus.maps = zeros(size(stimulus.rgbs));
    separatedrgbs = cell(size(stimulus.rgbs));
    for i = 1:length(stimulus.rgbs)
        rgb = stimulus.rgbs{i};
        map = find(sum(rgb ~= 0, 2));
        
        % This is weird, but agrees with the stimlisp output...
        if isempty(map), map = 0; end
        
        stimulus.maps(i) = map;
        
        separatedrgbs{i} = rgb(map,:);
    end
    stimulus.rgbs = separatedrgbs;
    
    stimulus.mapindices = unique(stimulus.maps);
end


% FIXME: Should really correct mapnyc as well, even though this is not
% currently used by allcones_plot
%
% Need to combine mapims and polys?
if length(stimulus.mapims) > 1
    stimulus.mapims = {combine_maps(stimulus.mapims)};
    stimulus.mapnycpoly = {cellfun(@(C)(C{1}), stimulus.mapnycpoly, 'UniformOutput', false)};
end