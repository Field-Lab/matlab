function stimstruct = parse_stim_rgbs(stimstruct, force)
% PARSE_STIM_RGBS   Categorize RGB stimuli for analysis convenience
%
% 2011 phli
% Redundant with some of Martin's functions
%

if nargin < 2, force = false; end
if isfield(stimstruct, 'urgbs') && ~force, return; end


if iscell(stimstruct.rgbs)
    [stimstruct.urgbs, ~, stimstruct.urgbi] = numcellunique(stimstruct.rgbs);
elseif isnumeric(stimstruct.rgbs)
    [stimstruct.urgbs, ~, stimstruct.urgbi] = numcellunique(num2cell(stimstruct.rgbs,2));
end

urgb.intensities = cell2mat(cellfun(@(a) (sum(a, 2)), stimstruct.urgbs, 'UniformOutput', false));
urgb.absolute_intensities = abs(urgb.intensities);
urgb.regions_stimulated  = urgb.absolute_intensities > 0;
urgb.nregions_stimulated = sum(urgb.regions_stimulated, 1);

urgb.blanks  = urgb.nregions_stimulated == 0;
urgb.singles = urgb.nregions_stimulated == 1;
urgb.doubles = urgb.nregions_stimulated == 2;

urgb.incr = cell2mat(cellfun(@(a)(all(a > 0, 2)), stimstruct.urgbs, 'UniformOutput', false));
urgb.decr = cell2mat(cellfun(@(a)(all(a < 0, 2)), stimstruct.urgbs, 'UniformOutput', false));

% Mixed polarity stimuli
urgb.uds = any(urgb.incr,1) & any(urgb.decr,1);

for i = 1:size(urgb.intensities, 1)
    uintensities{i} = unique(urgb.intensities(i,:));
end

stimstruct.urgb = urgb;
stimstruct.uintensities = uintensities;