function additivity_compound_plot(addstruct, rgc_indices, varargin)
% 2011-07 phli

opts = inputParser;
opts.addParamValue('az_pad_factor', 7);
opts.addParamValue('az_aspect_ratio', 1);
opts.addParamValue('mode', 'manhattan');
opts.KeepUnmatched = true;
opts.parse(varargin{:});
aropts = opts.Unmatched;
opts = opts.Results;

addrun = addstruct.addrun;
addstim = addstruct.addstim;
conerun = addstruct.conerun;
stabilityrun = addstruct.stablerun;

triggers = addrun.triggers;

if nargin < 2 || strcmp(rgc_indices, 'all')
    rgc_indices = 1:length(addrun.addrgcs);
end


for i = rgc_indices
    addrun_rgc_id  = addrun.addrgcs{i};
    if isempty(addrun_rgc_id), continue; end        
    additivity_raster_plot(addrun, addstim, addrun_rgc_id, 'subplot_width', 7, 'triggers', triggers, aropts);
    
    conerun_rgc_id = conerun.addrgcs(i);
    cones_1_and_2 = get_cones_by_weight_rank(conerun, conerun_rgc_id, [1 2]);

    rfax = subplot(2, 7, [4 5]);
    highlight_rgba = [1 0 0; 0 1 0];
    plot_voronoi_over_rf(conerun, conerun_rgc_id, 'highlight_cones', cones_1_and_2, 'highlight_rgba', highlight_rgba, opts);
    
    stableax = subplot(2, 7, [6 7]);
    overlay_cone_mosaics(conerun, stabilityrun, 'plot_borders', conerun, 'highlight_regions', cones_1_and_2, 'highlight_rgba', highlight_rgba, 'highlight_legend', {'1st' '2nd'});
    linkaxes([rfax, stableax])
end