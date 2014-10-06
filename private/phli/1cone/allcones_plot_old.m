function [figs rfweights rasterweights] = allcones_plot_old(datarun, conerun, varargin)
% ALLCONES_PLOT_OLD
% usage: allcones_plot(datarun, conerun, stimstruct, sorted_cone_indices, sorted_cone_weights, varargin)
%
% 2012-06, phli
%

stimstruct = regularize_allcones_stimulus(datarun.stimulus);
stimstruct = parse_stim_rgbs(stimstruct);

opts = inputParser();
opts.addParamValue('coneind', []);
opts.addParamValue('stablerun', []);
opts.addParamValue('mapweightindices', []);
opts.addParamValue('urgbs', stimstruct.urgbs);
opts.addParamValue('rfstimmapopts', {});
opts.addParamValue('stimmasks', cell(length(conerun.rgcs)));
opts.addParamValue('stimmaskcalcfunc', []);
opts.addParamValue('stimmasklocal', false);
opts.addParamValue('az_pad_factor', 6); % Used both for actually zooming in the RF plot and also for stimmasklocal calculation
opts.addParamValue('az_aspect_ratio', 1);
opts.addParamValue('LineWidth', 2);
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @jet);
opts.addParamValue('printpath', []);
opts.parse(varargin{:});
opts = opts.Results;

opts.rfstimmapopts = [opts.rfstimmapopts keepfields(opts, 'az_pad_factor', 'az_aspect_ratio', 'LineWidth')];

% Default local rectangle stim mask based on autozoom bounds
if opts.stimmasklocal, opts.stimmaskcalcfunc = @stimmasklocal; end

urgbs = opts.urgbs;
nrgbs = length(urgbs);
topheight = 1;
ploth = topheight+nrgbs;
nmaps = length(stimstruct.mapindices);

conestimsize = [conerun.stimulus.field_width conerun.stimulus.field_height];
stimmapsize = fliplr(size(stimstruct.mapims{1}));


figs = [];
rfweights = {};
rasterweights = {};
for i = 1:length(conerun.rgcs)
    rgc = datarun.rgcs{i};
    if isempty(rgc), continue; end
    clear rax lax sumhist sumhist2;
    polarity = conerun.stas.polarities{get_cell_indices(conerun, conerun.rgcs(i))};
    
    % Need to reload each time in case there was masking (see next step)
    stimmap = stimstruct.mapims{1};
    poly = stimstruct.mapnycpoly{1};
    
    % Mask the map (e.g. locally to exclude other parallel stimulated RGCs) if desired
    stimmask = opts.stimmasks{i};
    if ~isempty(opts.stimmaskcalcfunc)
        stimmask = opts.stimmaskcalcfunc(conerun, conerun.rgcs(i), stimmap, opts);
    end
    if ~isempty(stimmask)
        % Mask map, calculate new polygons
        stimmap = stimmap.*stimmask;
        poly = nyc2poly(map2manhattan(stimmap));
    end
    surviving_indices = sort(setdiff(unique(stimmap), 0)); % Unique should cover sort?    
    
    % Get predicted weights from integral of RF within stimulus region,
    rf_integral_weights = calc_stim_rf_weights(conerun, conerun.rgcs(i), stimmap);
    
    % Sort these so we can order raster plots by this calculated weight.
    % Must factor out the non-surviving indices from the sort.
    [sorted_rfiws, rfiw_indices] = sort(rf_integral_weights(surviving_indices), 'descend');
    rfiw_indices = surviving_indices(rfiw_indices);

    % Select specific stim regions after sorting by weight?
    if isempty(opts.mapweightindices)
        mapweightindices = 1:length(surviving_indices); 
    else
        mapweightindices = opts.mapweightindices;
        stimmap = selectmapindices(stimmap, rfiw_indices(mapweightindices));
        poly = nyc2poly(map2manhattan(stimmap));
        surviving_indices = sort(setdiff(unique(stimmap), 0)); % Unique should cover sort?    
    end
    nindices = length(surviving_indices);
    
    
    % Get colors to use, also calculate map index into colors
    if isempty(opts.colors), colors = opts.cmf(nindices);
    else colors = opts.colors; end
    mapcolors(surviving_indices,:) = colors;
    
    % Plot stim over RF from conerun first
    figs(i) = figure();
    sanesubplot(ploth, 3, {1:topheight 1});
    plot_stimmap_over_rf(conerun, conerun.rgcs(i), stimmap, 'poly', poly, 'colors', colors, opts.rfstimmapopts{:});
    cone_rf_bounds = axis();
    
    % Plot stability second if available
    if ~isempty(opts.stablerun)
        stablergc = opts.stablerun.rgcs{i};

        if ~isempty(stablergc)
            % Get bounds from cone rf plot and rescale
            stablestimsize = [opts.stablerun.stimulus.field_width opts.stablerun.stimulus.field_height];
            scale = stablestimsize ./ conestimsize;
            stable_bounds = cone_rf_bounds .* scale([1 1 2 2]);
            
            sanesubplot(ploth, 3, {1:topheight 2});
            plot_stimmap_over_rf(opts.stablerun, stablergc, stimmap, 'autozoom', false, 'bounds', stable_bounds, 'poly', poly, 'colors', colors, 'label_cones', false, opts.rfstimmapopts{:});
        end
    end
    
    % RF collapsed to greyscale within regions, for comparison to raster response greyscale
    sanesubplot(ploth, 3, {1 3});
    polycolors = rf_integral_weights(surviving_indices);
    polycolors = polycolors - min(polycolors);
    polycolors = polycolors ./ max(polycolors);
    if polarity == -1, polycolors = 1 - polycolors; end
    polycolors = repmat(polycolors(:), 1, 3);
    p = patchpolylines(poly, 'colors', colors, 'fillcolors', polycolors);
    set(gca, 'Color', 'w', 'Box', 'on', 'XTick', [], 'YTick', []);
    set(p, 'LineWidth', opts.LineWidth);
    axis ij;
    scale = stimmapsize ./ conestimsize;
    axis(cone_rf_bounds .* scale([1 1 2 2]));
    
    
    % FIXME: convert to URGB_RASTER_SUBPLOT format...
    surviving_coneids = recover_cones_stimulated(conerun, stimmap);
    for irgb = 1:length(urgbs);
        rgbs = stimstruct.rgbs;
        if ~iscell(rgbs) && isnumeric(rgbs)
            % Old pulse-sequence gets parsed into Nx3 RGB matrix instead of cell array of 1x3s
            rgbs = num2cell(rgbs, 2)';
        end
        rgbtrigs = cellfun(@(c)(isequal(c, urgbs{end-irgb+1})), rgbs);
        
        for imap = 1:nindices
            % Get mapindex according to rf stimulus integral weight order
            mapindex = rfiw_indices(mapweightindices(imap));
            iweight(irgb,imap) = rf_integral_weights(mapindex);
            
            sanesubplot(ploth, ceil(nindices*3/2), [topheight+irgb imap]);
            triggers = stimstruct.triggers(stimstruct.maps == stimstruct.mapindices(mapindex) & rgbtrigs);
            histcolor = mapcolors(mapindex,:);
            [~,ax,histx,hist{irgb,imap}] = rasterphli(datarun, rgc, triggers, 'hist', true, 'hist_color', histcolor, 'hist_bin', 0.05, 'plot', true, 'hist_line_width', 2, 'MarkerSize', 10, 'stop', mean(diff(stimstruct.triggers)));
            lax(irgb,imap) = ax(1);
            rax(irgb,imap) = ax(2);
            
            sumhist(irgb,imap) = sum(hist{irgb,imap});
            sumhist2(irgb,mapindex) = sumhist(irgb,imap);
            
            % Show headers above top row of rasters
            header = sprintf('r%.2f', sumhist(irgb,imap));
            if irgb == 1
                coneids = surviving_coneids{mapindex};
                header = sprintf('c %s\n%s', num2str(coneids), header);
            end
            title(header);
        end
        
        % Region shaded plot
        sanesubplot(ploth, 3, {topheight+irgb 3});
        polycolors = sumhist2(irgb, surviving_indices);
        polycolors = polycolors - min(polycolors);
        polycolors = polycolors ./ max(polycolors);
        if polarity == -1, polycolors = 1 - polycolors; end
        polycolors = repmat(polycolors(:), 1, 3);
        p = patchpolylines(poly, 'colors', colors, 'fillcolors', polycolors);
        set(gca, 'Color', 'w', 'Box', 'on', 'XTick', [], 'YTick', []);
        set(p, 'LineWidth', opts.LineWidth);
        axis ij;
        scale = stimmapsize ./ conestimsize;
        axis(cone_rf_bounds .* scale([1 1 2 2]));
        
%         % Mark-line plot
%         sanesubplot(ploth, 6, {topheight+irgb 6});
%         plot(iweight(irgb,:), sumhist(irgb,:), 'o-');
%         axis tight;
    end
    
    % Match y axes, remove superfluous ticks
    setaxesy(rax(:));
%     set(rax(:,1:end-1), 'YTick', []);

    % For now just take out all raster ticks below
%     % Remove superfluous #trials ticks if all rasters have same number of trials 
%     % Work by row
%     for r = 1:ploth
%         rlax = lax(r,:);
%         rlax = rlax(rlax > 0);
%         rlax_ylims = get(rlax, 'YLim');
%         if isequal(rlax_ylims{:})
%             set(rlax(2:end), 'YTick', []);
%         end
%     end
    
    % Remove superfluous time axes
    set(rax(1:end-1,:), 'XTick', []);
    lax = lax(lax > 0);
    set(lax(:), 'XTick', []);
    
    % Can be more subtle using above block, but right now that's slightly
    % messing up the plotting
    set(lax(:), 'YTick', []);
    
    % Save PDF
    if ischar(opts.printpath), print_plot(gcf, opts.printpath, datarun, conerun, i); end
    
    % Outputs
    rfweights{i}     = rf_integral_weights(surviving_indices);
    rasterweights{i} = sumhist2(:,         surviving_indices);
end

if nargout == 0, clear figs; end


function print_plot(fig, path, datarun, conerun, rgci)
parseconeprefix = parse_rrs_prefix(conerun);
parserunprefix = parse_rrs_prefix(datarun);
conergc = conerun.rgcs(rgci);
celltype = conerun.cell_types{find_cell_types(conerun, conergc)}.name;
name = sprintf('allcones_%s_%s_%s_f%02.fid%d',    ...
    parserunprefix.piece_fullname,          ...
    parserunprefix.run_full_name,           ...
    strrep(celltype, ' ', ''),              ...
    parseconeprefix.run_num,                ...
    conergc);

set(fig, 'PaperType', 'tabloid');
orient(fig, 'landscape');
print('-dpdf', fullfile(path, name));


% Assuming we're using default pad factor, etc., the mask should match the
% window shown on the stimmap over rf plot.
function mask = stimmasklocal(conerun, rgc, stimmap, opts)
mapsize = size(stimmap);
scale = mapsize ./ [conerun.stimulus.field_width conerun.stimulus.field_height];
bounds = autozoom_to_fit(conerun, rgc, opts.az_pad_factor, scale, opts.az_aspect_ratio);
rectmaskfunc = @(X,Y)(X > bounds(1) & X < bounds(2) & Y > bounds(3) & Y < bounds(4));
mask = funcmeshgrid(rectmaskfunc, 1:mapsize(2), 1:mapsize(1));


function map = selectmapindices(map, indices)
mask = map == indices(1);
for i = 2:length(indices)
    mask = mask | (map == indices(i));
end
map(~mask) = 0;