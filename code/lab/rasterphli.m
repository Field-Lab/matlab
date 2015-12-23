function [res,AX,histx,tt] = rasterphli(datarun, cell_specification, trigger, varargin)
% RASTERPHLI
% usage: [res,AX,histx,tt] = rasterphli(datarun,cell_specification,trigger,varargin)
%
% If HIST_BINS are given directly they must be evenly spaced.
%
%
% Derived from original by Martin
% phli 2011-07, does some silly (but helpful) things with histogram if there are no spikes in the window.
% phli 2012-02, reworked for modularity, may break pieces out, easier to run without plotting
%

p = inputParser;
p.addParamValue('plot', nargout == 0);
p.addParamValue('histax', 0);
p.addParamValue('rastax', 0);
p.addParamValue('start', 0);
p.addParamValue('stop', []);
p.addParamValue('stop_trigger_interval', 1);
p.addParamValue('color', [0 0 0]);
p.addParamValue('scale', 1);
p.addParamValue('MarkerStyle', '.');
p.addParamValue('MarkerSize', 10);
p.addParamValue('FontSize', 10);
p.addParamValue('axopts', {});
p.addParamValue('hist', nargout > 2);
p.addParamValue('hist_bin', .01);
p.addParamValue('hist_bins', []);
p.addParamValue('hist_color', [1 0 0]);
p.addParamValue('hist_line_width', 0.5);
p.addParamValue('LineStyle', '-');
p.addParamValue('YLim', []);
p.addParamValue('raster_ticks', length(trigger));
p.addParamValue('hist_ticks', false);
p.addParamValue('shadex', []);
p.addParamValue('shadecolor', 0.7*[1 1 1]);
p.addParamValue('plotopts', {});
p.addParamValue('pbaspect', []);
p.parse(varargin{:});
params = p.Results;


if ~isempty(params.hist_bins)
    params.start = params.hist_bins(1);
    params.stop  = params.hist_bins(end);
    
    hist_diff = diff(params.hist_bins);
    params.hist_bin = hist_diff(1);
end


indices = get_cell_indices(datarun, cell_specification);
spikes = cell2mat(datarun.spikes(indices));

if isempty(params.stop)
    params.stop = mean(diff(trigger)) * params.stop_trigger_interval;
end

res = [];
for i = 1:length(trigger),
    time_diffs = spikes - trigger(i);
    in_window = find(time_diffs >= params.start & time_diffs <= params.stop);
    res = [res; time_diffs(in_window), ones(size(in_window))*i];
end


if params.hist
    [histx, tt] = hist_rasters(res, trigger, params);
    params.histx = histx;
    params.tt = tt;
end

AX = [0 0];
if params.plot, AX = plot_rasters(res, trigger, params); end


if nargout == 0, clear res; end



function [histx, tt] = hist_rasters(res, trigger, params)
if isempty(params.hist_bins)
    params.hist_bins = params.start:params.hist_bin:params.stop;
end
histx = params.hist_bins;

if isempty(trigger)
    tt = NaN * ones(1, length(histx));
    return
end

if isempty(res)
    tt = zeros(1, length(histx));
    return
end

t = hist(res(:,1), histx);
tt = t / length(trigger) / params.hist_bin;

% FIXME: Ugly; can handle above instead?
% End bins are not full width
diffx = diff(histx);
half_bin_width = diffx(1) / 2;
first_bin_width = half_bin_width + histx(1)   - params.start;
last_bin_width  = half_bin_width - histx(end) + params.stop;
tt(1)   = tt(1)   * 2 * half_bin_width / first_bin_width;
tt(end) = tt(end) * 2 * half_bin_width / last_bin_width;



function AX = plot_rasters(res, trigger, params)
AX = [0 0];
if isempty(trigger), return; end

axopts = params.axopts;
axopts = [axopts 'FontSize', params.FontSize];

if params.rastax > 0, rastax = params.rastax;
else rastax = gca; end
AX(1) = rastax;

rastaxheight = (length(trigger) + 1) * params.scale;

% Useful for showing stimulus or boxcar template time range
if ~isempty(params.shadex)
    rectangle('Position', [params.shadex(1) 0 diff(params.shadex) rastaxheight], 'FaceColor', params.shadecolor, 'EdgeColor', 'none');
end

if ~isempty(res)
    plot(rastax, res(:,1),res(:,2)*params.scale, params.MarkerStyle, 'color', params.color, 'MarkerSize', params.MarkerSize);
end
axis(rastax, [params.start params.stop 0 rastaxheight]);
set(rastax, 'YTick', params.raster_ticks);
set(rastax, axopts{:});
if ~isempty(params.pbaspect), pbaspect(params.pbaspect); end


if params.hist
    if params.histax > 0, histax = params.histax;
    else histax = axes('Position', get(rastax, 'Position')); end
    AX(2) = histax;
    
    plot(histax, params.histx, params.tt, 'LineStyle', params.LineStyle, 'color', params.hist_color, 'LineWidth', params.hist_line_width);
    set(histax, 'YAxisLocation', 'right', 'XLim', get(rastax, 'XLim'), 'Color', 'none');
    set(histax, axopts{:});
    if ~isempty(params.pbaspect), pbaspect(params.pbaspect); end
    
    if ~isempty(params.YLim)
        set(histax, 'YLim', params.YLim); 
    else
        % At least don't plot with negatives...
        ylim = get(histax, 'YLim');
        set(histax, 'YLim', ylim(2)*[-0.02 1]);
    end
    
    if isnumeric(params.hist_ticks), set(histax, 'YTick', params.hist_ticks); end
end