function [rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds, varargin)

opts = inputParser();
opts.addParamValue('rasteropts', {});
opts.addParamValue('colors', []);
opts.addParamValue('cmf', []);
opts.parse(varargin{:});
opts = opts.Results;

if nargin < 6, bounds = []; end
if isempty(bounds), bounds = guess_bounds(conerun, rgcs); end

ncones = length(cones2plot);
colors = opts.colors;
if isempty(colors) && ~isempty(opts.cmf), colors = opts.cmf(ncones);       end
if isempty(colors),                       colors = repmat('k', ncones, 1); end

rasteropts = [opts.rasteropts 'hist_line_width' 2];
rfopts = {'fit', false, 'autozoom', false, 'scaled_up', 10};
triggers = rasterrun.triggers(1:2:end);
urgb = rasterrun.stimulus.urgb;
singleurgbs = find(urgb.singles);
nrgcs = length(rgcs);

outstructs = {};
figure();
for i = 1:nrgcs
    conerunid = rgcs(i);
    rasterrunid = map{get_cell_indices(conerun, conerunid)};
    
    % Load up URGB and HIST_COLOR with the proper trials and colors,
    % include space at top for RF plots.
    subplotxstart = (i-1) * 2;
    urgbs = cell(3,nrgcs*2);
    hist_color = urgbs;
    for j = 1:ncones
        cone = cones2plot(j);
        intensities = urgb.intensities(cone,singleurgbs);
        urgbs{j+2,1+subplotxstart} = singleurgbs(intensities ==  0.48*3);
        urgbs{j+2,2+subplotxstart} = singleurgbs(intensities == -0.48*3);
        hist_color{j+2,1+subplotxstart} = {colors(j,:)};
        hist_color{j+2,2+subplotxstart} = {colors(j,:)};
    end

    % Plot rasters for this RGC
    os = urgb_raster_subplot(rasterrun, rasterrunid, triggers, urgbs, 'color', {[1 1 1]./2}, 'hist_color', hist_color, rasteropts{:});
    outstructs = [outstructs os(3:end, subplotxstart+(1:2))];    
    
    % Plot RF for this RGC
    sanesubplot(2+ncones,nrgcs,{1:2 i});
    rfax(i) = plot_rf_stimmap(conerun, conerunid, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:}, 'scaled_up', 10);
    axis(bounds);
end


% Get bounds that include RF fit 2 sigma for all RGCs given
function bounds = guess_bounds(conerun, rgcs)
bounds = [Inf 0 Inf 0];
for cellid = rgcs
    cellnum = get_cell_indices(conerun, cellid); 
    
    fit = conerun.stas.fits{cellnum}; 
    mean = fit.mean; 
    sd = fit.sd;
    
    xmin = mean(1) - 2*sd(1); 
    xmax = mean(1) + 2*sd(1); 
    ymin = mean(2) - 2*sd(2); 
    ymax = mean(2) + 2*sd(2);
    
    bounds([1 3]) = min([bounds([1 3]); xmin ymin]); 
    bounds([2 4]) = max([bounds([2 4]); xmax ymax]);
end