function overlay_cone_mosaics(d1, d2, varargin)
% OVERLAY_CONE_MOSAICS  Show two cone mosaics overlaid, mostly for tracking movement
% usage: overlay_cone_mosaics(d1, d2, opts)
%
% opts
%   scale1              []
%   scale2              []
%   normalize_scale     true
%   plot_borders        []
%   label               [false false]
%
% 2011-07 phli
% See overlay_cone_mosaics_script for old code

opts = inputParser;
opts.addParamValue('scale1', []);
opts.addParamValue('scale2', []);
opts.addParamValue('bounds', []);
opts.addParamValue('normalize_scale', true);
opts.addParamValue('label', [false false]);
opts.addParamValue('clear', true);

opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;

if isempty(opts.scale1)
    if isfield(d1.cones, 'mosaic') && isfield(d1.cones.mosaic, 'voronoi_masks_scale') && ~isempty(d1.cones.mosaic.voronoi_masks_scale)
        opts.scale1 = d1.cones.mosaic.voronoi_masks_scale .* [1 1];
    else
        d1 = load_globals(d1);
        stim = stimulus_from_globals(d1.globals);
        opts.scale1 = [stim.stixel_width stim.stixel_height];
    end
end

if isempty(opts.scale2)
    if isfield(d1.cones, 'mosaic') && isfield(d2.cones.mosaic, 'voronoi_masks_scale') && ~isempty(d2.cones.mosaic.voronoi_masks_scale)
        opts.scale2 = d2.cones.mosaic.voronoi_masks_scale .* [1 1];
    else
        d2 = load_globals(d2);
        stim = stimulus_from_globals(d2.globals);
        opts.scale2 = [stim.stixel_width stim.stixel_height];
    end
end

if opts.normalize_scale
    opts.scale2 = opts.scale2 ./ opts.scale1;
    opts.scale1 = [1 1];
end


ax = gca;

plot_cone_mosaic(d1, 'fig_or_axes', ax, 'cone_size', 15, 'cone_colors', [0 0 0], 'scale', opts.scale1, 'bounds', opts.bounds, 'label', opts.label(1), 'label_color', [0 0 0], 'clear', opts.clear, unmatched)
plot_cone_mosaic(d2, 'fig_or_axes', ax, 'cone_size',  8, 'cone_colors', [1 1 1], 'scale', opts.scale2, 'bounds', opts.bounds, 'label', opts.label(2), 'label_color', [1 1 1], 'clear', false,      unmatched);