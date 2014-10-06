function ww_compound_plot(wwstruct, rgc_indices, varargin)
% WW_COMPOUND_PLOT
% usage: ww_compound_plot(wwstruct, rgc_indices, varargin)
%
% FIXME: Switch to using abstracted PLOT_STIMMAP_OVER_RF
%
% 2011 phli
%

opts = inputParser;
opts.addParamValue('cmf', @jet);
opts.addParamValue('colors', [1 0 0; 0 0.5 0; 0 0 1; 1 0 1]);
opts.addParamValue('az_pad_factor', 4);
opts.addParamValue('az_aspect_ratio', 1);
opts.addParamValue('plot_voronoi', false);
opts.addParamValue('scaled_up', []);
opts.addParamValue('printpath', []);
opts.addParamValue('closefigs', false);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
wwropts = opts.Unmatched;
opts = opts.Results;

wwrun = wwstruct.wwrun;
wwstim = wwrun.stimulus;
conerun = wwstruct.conerun;

stabilityrun = [];
if isfield(wwstruct, 'stablerun')
    stabilityrun = wwstruct.stablerun;
end

if nargin < 2 || strcmp(rgc_indices, 'all')
    rgc_indices = 1:length(wwrun.wwrgcs);
end

% Default opts for printing
if isempty(opts.scaled_up), 
    if ischar(opts.printpath)
    	opts.scaled_up = 10;
    else
        opts.scaled_up = 1; 
    end
end

if ~isfield(wwropts, 'MarkerSize'), wwropts.MarkerSize = 15; end
if ~isfield(wwropts, 'FontSize'),   wwropts.FontSize   = 14; end

wwropts.hist_colors = opts.colors;
wwropts.cmf = opts.cmf;
wwropts.subplot_add_width = 2;
wwropts.subplot_add_height = 0;

for i = 1:length(rgc_indices)
    rgc = wwrun.wwrgcs(rgc_indices(i));
    if isempty(rgc{1}), continue; end

    conergc = conerun.wwrgcs(rgc_indices(i));
    if iscell(conergc), conergc = conergc{1}; end
    polarity = get_polarities(conerun, conergc);
    
    ww_raster_plot(wwrun, wwstim, rgc, polarity, wwropts);
    
    if opts.plot_voronoi && isfield(conerun, 'cones')
        ww_rf_voronoi_stim_plots(conerun, stabilityrun, wwstim, rgc_indices(i), opts);
    else
        ww_rf_stim_plots(conerun, stabilityrun, wwstim, rgc_indices(i), opts);
    end
    
    if ischar(opts.printpath)
        mkdir(opts.printpath);
        print_plot(gcf, opts.printpath, wwrun, conerun, rgc_indices(i));
    end
    
    if opts.closefigs, close; end
end


function print_plot(fig, path, wwrun, conerun, rgci)
name = build_1cone_filename(conerun, wwrun, conerun.wwrgcs(rgci));
name = ['ud_' name];
set(fig, 'PaperType', 'tabloid');
orient(fig, 'landscape');
print('-dpdf', fullfile(path, name));



function ax1 = ww_rf_stim_plots(conerun, stabilityrun, wwstim, rgc_index, opts)
rgc = conerun.wwrgcs(rgc_index);
if iscell(rgc), rgc = rgc{1}; end

ax1 = sanesubplot(4, 6, {1:2 5:6});
optnames = {'az_pad_factor', 'az_aspect_ratio', 'scaled_up', 'colors', 'cmf'};
plotopts = keepfields(opts, optnames{:}, '-force');
plotopts.fit = false;
plotopts.patchopts = {'LineWidth', 2};
plot_rf_stimmap(conerun, rgc, wwstim.mapnycpoly{1}, plotopts);

rgc2 = [];
if isfield(stabilityrun, 'wwrgcs'), rgc2 = stabilityrun.wwrgcs{rgc_index}; end
if ~isempty(rgc2) % No match
    ax2 = sanesubplot(4, 6, {3:4 5:6});
    plot_rf_stimmap(stabilityrun, rgc2, wwstim.mapnycpoly{1}, plotopts);
end



function ax1 = ww_rf_voronoi_stim_plots(conerun, stabilityrun, wwstim, rgc_index, opts)
rgc = conerun.wwrgcs(rgc_index);
if iscell(rgc), rgc = rgc{1}; end

bounds = autozoom_to_fit(conerun, rgc, opts.az_pad_factor, 1, opts.az_aspect_ratio);
ax1 = sanesubplot(4, 6, {1:2 5:6});
plot_voronoi_over_rf(conerun, rgc, 'bounds', bounds, 'scaled_up', opts.scaled_up);
scalemap = 1 / conerun.cones.mosaic.voronoi_masks_scale;
hold on;
patchmanhattan(wwstim.mapnyc{1}, 'scale', [scalemap scalemap], 'patchopts', {'LineWidth', 2}, 'colors', opts.colors, 'cmf', opts.cmf);
autozoom_to_fit(conerun, rgc, opts.az_pad_factor, 1, opts.az_aspect_ratio);

rgc2 = [];
if isfield(stabilityrun, 'wwrgcs'), rgc2 = stabilityrun.wwrgcs{rgc_index}; end
ax2 = sanesubplot(4, 6, {3:4 5:6});
if isempty(rgc2) % No match
    bounds = autozoom_to_fit(conerun, rgc, opts.az_pad_factor * 2, 1, opts.az_aspect_ratio); % A little bigger
    scale2 = stabilityrun.cones.mosaic.voronoi_masks_scale ./ conerun.cones.mosaic.voronoi_masks_scale;
    overlay_cone_mosaics(conerun, stabilityrun, 'scale1', [1 1], 'scale2', [scale2 scale2], 'bounds', stabilitybounds);
else
    bounds = autozoom_to_fit(stabilityrun, rgc2, opts.az_pad_factor, 1, opts.az_aspect_ratio);
    plot_voronoi_over_rf(stabilityrun, rgc2, 'bounds', bounds, 'scaled_up', opts.scaled_up);
    scalemap = 1 / stabilityrun.cones.mosaic.voronoi_masks_scale;
end
hold on;
patchmanhattan(wwstim.mapnyc{1}, 'scale', [scalemap scalemap], 'patchopts', {'LineWidth', 2},  'colors', opts.colors, 'cmf', opts.cmf);
axis(bounds);