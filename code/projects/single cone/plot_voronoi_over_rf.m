function plot_voronoi_over_rf(datarun, rgc, varargin)
% 2011-07 phli

opts = inputParser;
opts.addParamValue('voronoirun', datarun);
opts.addParamValue('title', true);
opts.addParamValue('color_transform', []);
opts.addParamValue('scaled_up', 1);
opts.addParamValue('autozoom', true);
opts.addParamValue('az_pad_factor', 5);
opts.addParamValue('az_aspect_ratio', []);
opts.addParamValue('bounds', []);
opts.addParamValue('mode', 'polygons');
opts.addParamValue('plot_cone_centers', false);
opts.addParamValue('border_color', 'k');
opts.addParamValue('label_cones', false);
opts.addParamValue('label_color', 'y');
opts.addParamValue('highlight_cones', []);
opts.addParamValue('highlight_ranked_cones', []);
opts.addParamValue('highlight_rgba', []);
opts.addParamValue('cmf', @jet);
opts.addParamValue('highlight_legend', {});
opts.parse(varargin{:});
opts = opts.Results;


% Plot RF
plot_rf(datarun, rgc, 'fit', false, 'ticks', false, 'scale', opts.scaled_up, 'title', opts.title, 'color_transform', opts.color_transform);
hold on;

if opts.autozoom
    opts.bounds = autozoom_to_fit(datarun, rgc, opts.az_pad_factor, [1 1], opts.az_aspect_ratio);
end

pvopts = keepfields(opts, 'mode', 'bounds', 'highlight_rgba', 'highlight_legend', 'highlight_ranked_cones');
pvopts.plot_centers = opts.plot_cone_centers;
pvopts.border_color = opts.border_color;
pvopts.label_cones = opts.label_cones;
pvopts.label_color = opts.label_color;
pvopts.highlight_regions = opts.highlight_cones;
pvopts.cone_rank_rgc = rgc;
plot_voronoi(opts.voronoirun, pvopts);