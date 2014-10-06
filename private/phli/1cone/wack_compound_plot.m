function wack_compound_plot(wackstruct, rgc_indices, varargin)

opts = inputParser;
opts.addParamValue('az_pad_factor', 3);
opts.addParamValue('az_aspect_ratio', 1);
opts.addParamValue('mode', 'manhattan');
opts.KeepUnmatched = true;
opts.parse(varargin{:});
wropts = opts.Unmatched;
opts = opts.Results;

wackrun = wackstruct.wackrun;
cells = wackrun.wackrgcs;
wackstim = wackstruct.wackstim;
conerun = wackstruct.conerun;
stabilityrun = wackstruct.stablerun;

if nargin < 2 || strcmp(rgc_indices, 'all')
    rgc_indices = 1:length(cells);
end

if isfield(wackrun, 'rgc_regions')
    cellregions = wackrun.rgc_regions;
else
    cellregions = build_wack_regions(wackstruct);
end

wropts.subplot_width = 6;
wack_raster_plot(wackrun, wackstim, cells(rgc_indices), cellregions(rgc_indices), wropts);


for i = 1:length(rgc_indices)
    rgc_index = rgc_indices(i);
    conerun_rgc_id = conerun.conergcs(rgc_index);
    cones_1_and_2 = get_cones_by_weight_rank(conerun, conerun_rgc_id, [1 2]);

    rfax = subplot(length(rgc_indices), wropts.subplot_width, sub2ind([wropts.subplot_width length(rgc_indices)], 5, i));
    highlight_rgba = [1 0 0; 0 1 0];
    plot_voronoi_over_rf(conerun, conerun_rgc_id, 'highlight_cones', cones_1_and_2, 'highlight_rgba', highlight_rgba, opts);
    
    stableax = subplot(length(rgc_indices), wropts.subplot_width, sub2ind([wropts.subplot_width length(rgc_indices)], 6, i));
    overlay_cone_mosaics(conerun, stabilityrun, 'plot_borders', conerun, 'highlight_regions', cones_1_and_2, 'highlight_rgba', highlight_rgba, 'highlight_legend', {'1st' '2nd'});
    linkaxes([rfax, stableax])
end


% Assumes that the regions were numbered as RGC1_cone1, RGC1_cone2,
% RGC2_cone1, etc.  May have to make this more flexible to handle
% older/future schemes...
function cellregions = build_wack_regions(wackstruct)

numregions = size(wackstruct.wackstim.intensities, 1);
numrgcs = length(wackstruct.conerun.conergcs);

if numregions ~= 2*numrgcs
    error('Wack cannot infer voronoi regions for each RGC; # map indices is not 2*#RGCs');
end

% [1, 3, 5, 7 ...]
firstregions = 1:2:numregions;

% Converts to {1:2 3:4 5:6 7:8 ...}
cellregions = cellfun(@(e) (e:e+1), num2cell(firstregions), 'UniformOutput', false);