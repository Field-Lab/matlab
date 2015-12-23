function [crs rasterhistsy] = single_map_contrast_response(datarun, triggers, rgcs, template, urgbs, varargin)
% SINGLE_MAP_CONTRAST_RESPONSE  Template based contrast response calc
% usage: [crs rasterhistsy] = single_map_contrast_response(datarun, triggers, rgcs, template, urgbs, varargin)
%
% Inefficient implementation, but getting the job done.  This is for
% stimuli with a single map file with multiple indices, so it assumes that
% stimrgbs is the only stimulus information that we need to operate on to
% distinguish different trials.
%
% 2012-02 phli
%

crs = nan(length(rgcs), length(urgbs));
rasterhistsy = cell(length(rgcs), length(urgbs));
for j = 1:length(urgbs)
    urgb = urgbs(j);
    
    % URGBS can be either a cell array of stimulus values to match, or just
    % a numeric array of indices into DATARUN.STIMULUS.URGBI.
    if iscell(urgb)
        % Find the STIMULUS.URGB that matches the given stimulus values
        urgb = urgb{1};
        urgbmatching = cellfun(@(a)(isequal(a,urgb)), datarun.stimulus.urgbs);
        urgb = find(urgbmatching);
    end
    selected_trials = datarun.stimulus.urgbi{urgb};

    % Take only selected trials that were actually run (based on number of triggers)
    selected_trials = selected_trials(selected_trials < length(triggers));
    selected_triggers = triggers(selected_trials);
    if isempty(selected_triggers), continue; end
    
    for k = 1:length(rgcs)
        if iscell(rgcs), rgc = rgcs{k};
        else rgc = rgcs(k); end
        if isempty(rgc), continue; end
        
        [~, ~, ~, rasterhisty] = rasterphli(datarun, rgc, selected_triggers, 'hist_bins', template.x, varargin{:});
        
        crs(k,j) = template.y(:)' * rasterhisty(:);
        rasterhistsy{k,j} = rasterhisty;
    end
end