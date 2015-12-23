function [rasterrun, conerun, stabilityrun, crs, crsx, ps, resnorms, residuals] = cr_compound_plot(rasterrun, conerun, triggers, varargin)
% CR_COMPOUND_PLOT 
%
% Currently does not handle using different C/R templates for up versus
% down stimuli.
%
% 2012-10, phli
%

opts = inputParser();
opts.addParamValue('rgcindices', 'all');
opts.addParamValue('simple', false); % Hack to handle situations with only 1 direction of contrasts (e.g. in allcones analyses)
opts.addParamValue('usedmapindices', {});  % Another special case for handling allcones analyses
opts.addParamValue('blankmapindices', {}); % Another special case for handling allcones analyses
opts.addParamValue('template', []);
opts.addParamValue('baselines', []);
opts.addParamValue('boxstart', []);
opts.addParamValue('boxend', 0.25);
opts.addParamValue('shadex', []);
opts.addParamValue('cmf', @jet);
opts.addParamValue('colors', []);
opts.addParamValue('stabilityrun', []);
opts.addParamValue('printpath', []);
opts.addParamValue('closefigafterprint', true);
opts.addParamValue('scaled_up', []);
opts.addParamValue('surround', false);
opts.addParamValue('title', false);
opts.addParamValue('marksthresh', 3.5);
opts.addParamValue('YLim', true);
opts.parse(varargin{:});
opts = opts.Results;

if strcmp(opts.rgcindices, 'all'), opts.rgcindices = 1:length(rasterrun.rgcs); end
if isempty(opts.boxstart)
    if opts.surround, opts.boxstart = 0.075;
    else              opts.boxstart = 0.05; end
end

if isempty(opts.scaled_up)
    if ischar(opts.printpath), opts.scaled_up = 10;
    else                       opts.scaled_up = 1; end
end

% Defaults to boxcar template
if isempty(opts.template)
    template.x = -0.1:0.02:0.5;
    template.y = zeros(size(template.x));
    template.y(template.x > opts.boxstart & template.x < opts.boxend) = 1;
    template.y = template.y ./ sum(template.y);
    opts.shadex = [opts.boxstart opts.boxend];
else
    template = opts.template;
end

% Load if needed
rasterrun = read_stim_lisp_output(rasterrun);
rasterrun.stimulus = parse_stim_rgbs(rasterrun.stimulus);
rasterrun.stimulus = parse_cr_rgbs(rasterrun.stimulus);
stimulus = rasterrun.stimulus;
if ~isfield(conerun, 'stas')
    conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
    conerun = set_polarities(conerun);
end
if ~isempty(opts.stabilityrun) && ~isfield(opts.stabilityrun, 'stas')
    opts.stabilityrun = load_sta(opts.stabilityrun, struct('load_sta', [], 'guess_stimulus', false));
    opts.stabilityrun = set_polarities(opts.stabilityrun);
end

% Precalc the whole block of C/R values, rearrange into convenient dimensionality
[crs rasterhistsy] = single_map_contrast_response(rasterrun, triggers, rasterrun.rgcs, template, cell2mat(stimulus.cr.single_cones'));
crs = permute(reshape(crs, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
rasterhistsy = permute(reshape(rasterhistsy, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
crsx = cell2mat(rasterrun.stimulus.cr.single_cone_intensities);
ncontrasts = size(crs,2); % Assume same for all

% Setup plotting params/styles
rastercoltitles = cellfun(@num2str, num2cell(stimulus.cr.single_cone_intensities{1} ./ 3), 'UniformOutput', false); % Assume same for all

% Plot paneled figures for each RGC
for i = opts.rgcindices
    baseline = [];
    if ~isempty(opts.baselines) && length(opts.baselines) >= i, baseline = opts.baselines(i); end

    
    % Hacks for allcones analyses
    usedmapindices = 1:length(stimulus.cr.single_cones);
    if ~isempty(opts.usedmapindices) && length(opts.usedmapindices) >= i, usedmapindices = opts.usedmapindices{i}; end    
    if isempty(opts.colors), opts.colors = opts.cmf(length(usedmapindices)); end
    histcolors = mat2cell(opts.colors, ones(size(opts.colors,1), 1));
    
    % Special case for allcones analyses missing blanks; use 
    % unused map indices as the blanks.  Just to be safe, use the
    % lowest contrast stimulus and average
    blankmapindices = [];
    if ~isempty(opts.blankmapindices) && length(opts.blankmapindices) >= i, blankmapindices = opts.blankmapindices{i}; end
    if isempty(baseline) && ~isempty(blankmapindices)
        blankcrsx = crsx(blankmapindices,:);
        blankcrs  = crs(blankmapindices, :);
        lowc = abs(blankcrsx) == min(abs(blankcrsx(:)));
        baseline = mean(blankcrs(lowc));
    end
    
    
    rgc = rasterrun.rgcs{i};
    if isempty(rgc), continue; end
    figure(rgc);
    
    % Plot C/R rasters (assume same number in each row)
    rasterploturgbs = num2cell(cell2mat(stimulus.cr.single_cones(usedmapindices)));
    urgb_raster_subplot(rasterrun, rgc, triggers, rasterploturgbs, 'add_width', 2, 'shadex', opts.shadex, 'titles', rastercoltitles, 'hist_color', histcolors, 'YLim', opts.YLim);
    
    % How many rows in last col?
    numrows = 2;
    if ~isempty(opts.stabilityrun), numrows = 3; end
    
    % Get conerun RGC, robust to either numeric or cell array
    conergc = conerun.rgcs(i);
    if iscell(conergc), conergc = conergc{1}; end
    
    % Stability plots
    sanesubplot(numrows, ncontrasts+2, {1 (1:2)+ncontrasts});
    rfopts = keepfields(opts, 'colors', 'scaled_up');
    if opts.surround, rfopts.az_pad_factor = 15; end
    conerun = setsimplemarks(conerun, conergc, opts.marksthresh);
    plot_rf_stimmap(conerun, conergc, stimulus.mapnycpoly{1}(usedmapindices), 'fit', false, 'pdtitle', false, rfopts);
    if ~isempty(opts.stabilityrun) && ~isempty(opts.stabilityrun.rgcs{i})
        sanesubplot(numrows, ncontrasts+2, {2 (1:2)+ncontrasts});
        opts.stabilityrun = setsimplemarks(opts.stabilityrun, opts.stabilityrun.rgcs{i}, opts.marksthresh);
        plot_rf_stimmap(opts.stabilityrun, opts.stabilityrun.rgcs{i}, stimulus.mapnycpoly{1}(usedmapindices), 'fit', false, 'pdtitle', false, rfopts);
    end
    
    % Plot xscaled C/R
    sanesubplot(numrows, ncontrasts+2, {numrows (1:2)+ncontrasts});
    if opts.simple
        [ps{i} resnorms(i) residuals{i}] = normcdfxscalesimple(crs(usedmapindices,:,i), crsx(usedmapindices,:), 'plot', true, 'baseline', baseline, 'subtract_baseline', true, 'colors', opts.colors, 'title', opts.title);
    else
        polarity = get_polarities(conerun, conergc);
        if opts.surround, polarity = polarity * -1; end
        [ps{i} resnorms(i) residuals{i}] = normcdfxscale(crs(:,:,i), crsx, 'plot', true, 'baseline', baseline, 'polarity', polarity, 'rectfit', false, 'subtract_baseline', true, 'colors', opts.colors, 'title', opts.title);
    end
    axis square
    
    % Save PDF?
    if ischar(opts.printpath)
        print_plot(gcf, opts.printpath, rasterrun, conerun, i);
        if opts.closefigafterprint, close(gcf); end
    end
end


if nargout == 0, clear conerun; end
if nargout > 2, stabilityrun = opts.stabilityrun; end


function print_plot(fig, path, rasterrun, conerun, rgci)
figure(fig);
drawnow();

parseconeprefix = parse_rrs_prefix(conerun);
parserunprefix = parse_rrs_prefix(rasterrun);
conergc = conerun.rgcs(rgci);
if iscell(conergc), conergc = conergc{1}; end
celltype = conerun.cell_types{find_cell_types(conerun, conergc)}.name;
name = sprintf('cr_%s_%s_%s_f%02.fid%d',    ...
    parserunprefix.piece_fullname,          ...
    parserunprefix.run_full_name,           ...
    strrep(celltype, ' ', ''),              ...
    parseconeprefix.run_num,                ...
    conergc);

set(fig, 'PaperType', 'tabloid');
orient(fig, 'landscape');
print('-dpdf', fullfile(path, name));