function h = plot_stimmap_over_rf(rfrun, rgc, mapim, varargin)
% PLOT_STIMMAP_OVER_RF
% usage: h = plot_stimmap_over_rf(rfrun, rgc, map, opts)
%
% 2012-07, phli
%


opts = inputParser;

% For plot_rf
opts.addParamValue('title', true);
opts.addParamValue('scaled_up', 1);

% For plotting stimuli
opts.addParamValue('mode', 'polygons');
opts.addParamValue('colors', []);
opts.addParamValue('LineWidth', 2);
opts.addParamValue('cmf', @jet)
opts.addParamValue('stimopts', {});

% Can pass these if already precalculated
opts.addParamValue('nyc', []);
opts.addParamValue('poly', []);

% Label cones
opts.addParamValue('load_cones_index', []);
opts.addParamValue('label_cones', true);
opts.addParamValue('label_color_match_stimcolor', true);
opts.addParamValue('label_unstimulated_cones', false);
opts.addParamValue('label_default_color', 'y');

% Bounds
opts.addParamValue('bounds', [0 rfrun.stimulus.field_height 0 rfrun.stimulus.field_width]);
opts.addParamValue('autozoom', true);
opts.addParamValue('az_pad_factor', 5);
opts.addParamValue('az_aspect_ratio', []);

opts.parse(varargin{:});
opts = opts.Results;


if ~isfield(rfrun, 'stimulus')
    rfrun = load_sta(rfrun, 'load_sta', rgc);
end
rfrun = set_polarities(rfrun);
rfrun = get_sta_fits_from_vision(rfrun, rgc);



% Plot RF
plot_rf(rfrun, rgc, 'fit', false, 'ticks', false, 'scale', opts.scaled_up, 'title', opts.title);
hold on;

if opts.autozoom, 
    opts.bounds = autozoom_to_fit(rfrun, rgc, opts.az_pad_factor, [1 1], opts.az_aspect_ratio);
end



% Calc number of stimuli
stimindices = setdiff(unique(mapim), 0);
if isempty(stimindices) 
    h = [];
    axis(opts.bounds);
    return
end
nstims = length(stimindices);

% Calculate colors for stim
if isempty(opts.colors), opts.colors = opts.cmf(nstims); end
opts.stimopts = [opts.stimopts 'colors' opts.colors];

% Calculate scaling needed for stim
scalemap = [1 1] ./ [rfrun.stimulus.stixel_width rfrun.stimulus.stixel_height];
opts.stimopts = [opts.stimopts 'scale' scalemap];

% Pick stim plot mode
switch opts.mode
    case 'polygons'
        plotfunc = @patchpolylines;
        if isempty(opts.poly)
            if isempty(opts.nyc), opts.nyc = map2manhattan(mapim); end
            opts.poly = nyc2poly(opts.nyc);
        end
        plotdata = opts.poly;
    case 'segs'
        plotfunc = @patchmanhattan;
        if isempty(opts.nyc), opts.nyc = map2manhattan(mapim); end
        plotdata = opts.nyc;
    case 'fills'
        plotfunc = @fillmap;
        plotdata = mapim;
end

% Plot stim
h = plotfunc(plotdata, opts.stimopts{:});
set(h, 'LineWidth', opts.LineWidth);
if nargout < 1, clear h; end



% Plot cone labels?
if opts.label_cones
    if ~isfield(rfrun, 'cones'), rfrun = load_cones(rfrun, opts.load_cones_index); end
    centers = rfrun.cones.centers;
    ncones = size(centers,1);
    cids = 1:ncones;
    
    % Determine colors for labels
    colors = repmat({opts.label_default_color}, [ncones, 1]);
    recovered_cones = recover_cones_stimulated(rfrun, mapim);
    if opts.label_color_match_stimcolor
        for i = 1:length(opts.colors)
            if stimindices(i) > length(recovered_cones), continue; end
            colors(recovered_cones{stimindices(i)}) = {opts.colors(i,:)};
        end
    end
    
    if ~opts.label_unstimulated_cones
        recovered_cones = cell2mat(recovered_cones(:));
        centers = centers(recovered_cones,:);
        cids    = cids(recovered_cones);
        colors  = colors(recovered_cones);
    end
    
    inxbounds = centers(:,1) > opts.bounds(1) & centers(:,1) < opts.bounds(2);
    inybounds = centers(:,2) > opts.bounds(3) & centers(:,2) < opts.bounds(4);
    inbounds = inxbounds & inybounds;
    
    for i = 1:size(centers,1)
        if ~inbounds(i), continue; end
        t = text(centers(i,1), centers(i,2), num2str(cids(i)));
        set(t, 'Color', colors{i});
    end
end



% Zoom
axis(opts.bounds);