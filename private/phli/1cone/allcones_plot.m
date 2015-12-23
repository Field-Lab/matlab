function [ax rfweights crweights f gof p crs crsx] = allcones_plot(conerun, rasterrun, stimmap, padfactor, varargin)

opts = inputParser();
opts.addParamValue('ar', 1);
opts.addParamValue('stablerun', struct([]));
opts.addParamValue('rfopts',     {});
opts.addParamValue('xscaleopts', {});
opts.addParamValue('marksthresh', 3.5);
opts.parse(varargin{:});
opts = opts.Results;

regions = setdiff(unique(stimmap), 0);

mapsize = size(stimmap);
scale = mapsize ./ [conerun.stimulus.field_width conerun.stimulus.field_height];
bounds = autozoom_to_fit(conerun, conerun.rgcs, padfactor, scale, opts.ar);
localpoly = nyc2poly(map2manhattan(stimmap));

% Plot RF
ax(1,1) = subplot(3,2,1);
conerun = setsimplemarks(conerun, conerun.rgcs, opts.marksthresh);
plot_rf_stimmap(conerun, conerun.rgcs, localpoly, 'fit', false, 'az_pad_factor', padfactor, 'az_aspect_ratio', opts.ar, opts.rfopts{:});

% Plot stability?
if ~isempty(opts.stablerun) && ~isempty(opts.stablerun.rgcs{1})
    ax(1,2) = subplot(3,2,2);
    opts.stablerun = setsimplemarks(opts.stablerun, opts.stablerun.rgcs{1}, opts.marksthresh);
    plot_rf_stimmap(opts.stablerun, opts.stablerun.rgcs{1}, localpoly, 'fit', false, 'az_pad_factor', padfactor, 'az_aspect_ratio', opts.ar, opts.rfopts{:});
end

% RF weights
rfweights = calc_stim_rf_weights(conerun, conerun.rgcs, stimmap);
rfweights = rfweights(regions);

% CR fits
[crs crsx rasterhistx rasterhisty] = calc_allcones_cr(rasterrun, rasterrun.rgcs{1}, rasterrun.triggers(1:2:end));
ax(2,1) = subplot(3,2,3);
[p resnorm residual] = normcdfxscalesimple(crs(regions,:), crsx(regions,:), 'plot', true, 'title', false, opts.xscaleopts{:});
axis square
crweights = p(1:end-2)';


% Calculate normalization, check
rfweights = rfweights./max(rfweights);
[f gof] = fit(crweights, rfweights, fittype({'x'}));
ax(2,2) = subplot(3,2,4);
plot(rfweights, f.a.*crweights, '.k'); hold on; plot([0 1], [0 1], 'k');
axis square

% Plot allcones!
maxweight = max([max(rfweights) max(f.a*crweights)]);
minweight = min([min(rfweights) min(f.a*crweights)]);
ax(3,1) = subplot(3,2,5);
patchpolylines(localpoly, 'colors', repmat('k', length(regions), 1), 'fillcolors', repmat((maxweight-rfweights)./(maxweight-minweight), 1, 3)); 
axis(bounds); axis ij
ax(3,2) = subplot(3,2,6);
patchpolylines(localpoly, 'colors', repmat('k', length(regions), 1), 'fillcolors', repmat((maxweight-f.a*crweights)./(maxweight-minweight), 1, 3));
axis(bounds); axis ij