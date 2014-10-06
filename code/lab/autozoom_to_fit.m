function bounds = autozoom_to_fit(datarun, cell_id, pad_factor, scale, aspect_ratio)
% AUTOZOOM_TO_FIT   Zoom in on plot based on RF fit
% usage: bounds = autozoom_to_fit(datarun, cell_id, pad_factor, scale, aspect_ratio)
%
% Output is [xstart xend ystart yend]
%
% 2011-07 phli


if nargin < 3
    pad_factor = 5;
end

if nargin < 4
    scale = [1 1];
end
if length(scale) == 1
    scale = [scale scale];
end

if nargin < 5
    aspect_ratio = 1;
end


cell_num = get_cell_indices(datarun, cell_id);
fit = datarun.stas.fits{cell_num};

xdiff = fit.sd(1) / 2 * pad_factor;
ydiff = fit.sd(2) / 2 * pad_factor;

if ~isempty(aspect_ratio)
    ydiff = max(ydiff, xdiff/aspect_ratio);
    xdiff = max(xdiff, ydiff*aspect_ratio);
end

xstart = fit.mean(1) - xdiff;
xend   = fit.mean(1) + xdiff;
ystart = fit.mean(2) - ydiff;
yend   = fit.mean(2) + ydiff;
bounds = [scale(1) scale(1) scale(2) scale(2)].*[xstart xend ystart yend];

if nargout < 1
    if ~any(isnan(bounds)), axis(bounds); end
    clear bounds;
end