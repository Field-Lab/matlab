function stimstruct = parse_cr_rgbs(stimstruct, blanks)
% PARSE_CR_RGBS   Categorize contrast response stimuli for analysis convenience
%
% 2012 phli
%

if nargin < 2, blanks = true; end

% Do we have what we need to work with?
if ~isfield(stimstruct, 'urgbs'), stimstruct = parse_stim_rgbs(stimstruct); end

% Get the single cone stimulated trials for each cone, sorted by contrast
cr.single_cones            = cell(stimstruct.numcones,1);
cr.single_cone_intensities = cell(stimstruct.numcones,1);
for i = 1:stimstruct.numcones
    stims = stimstruct.urgb.singles & stimstruct.urgb.regions_stimulated(i,:);
    if blanks, stims = stims | stimstruct.urgb.blanks; end
    stims = find(stims);
    
    intensities = stimstruct.urgb.intensities(i,stims);
    [~,ind] = sort(intensities);
    cr.single_cones{i}            = stims(ind);
    cr.single_cone_intensities{i} = intensities(ind);
end

% Get the double cone stimulated trials for each pair, sorted by contrast
% (In current data through 2012-10 these are not that useful as we didn't collect many contasts)
crdoubles = stimstruct.urgb.doubles & ~stimstruct.urgb.uds;
cr.double_cones = {};
cr.double_cone_intensities = {};
for i = 1:stimstruct.numcones
    for j = (i+1):stimstruct.numcones
        stims = crdoubles & stimstruct.urgb.regions_stimulated(i,:) & stimstruct.urgb.regions_stimulated(j,:);
        if blanks, stims = stims | stimstruct.urgb.blanks; end
        stims = find(stims);
        
        intensities = sum(stimstruct.urgb.intensities([i j],stims));
        [~,ind] = sort(intensities);
        cr.double_cones{end+1,1}            = stims(ind);
        cr.double_cone_intensities{end+1,1} = intensities(ind);
    end
end

stimstruct.cr = cr;