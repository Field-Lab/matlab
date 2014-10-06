function hout = sliderplot_ei_contours(contours, varargin)

opts = inputParser();
opts.addParamValue('xlims', [-360 360]);
opts.addParamValue('ylims', [-390 390]);
opts.addParamValue('zlims', []);
opts.parse(varargin{:});
opts = opts.Results;

f = figure();

if isempty(opts.zlims)
    minelevation = min(cell2mat(collect(contours, @(c) min([c.elevation]))));
    maxelevation = max(cell2mat(collect(contours, @(c) max([c.elevation]))));
    opts.zlims = [minelevation maxelevation];
end
setappdata(f, 'xlims', opts.xlims);
setappdata(f, 'ylims', opts.ylims);
setappdata(f, 'zlims', opts.zlims);


gspopts.data = contours;
gspopts.handle = f;
gspopts.data_length_func = @length;
gspopts.get_imag_func    = @get_ei_contour_frame;
gspopts.plot_func        = @plot_func;
fig = gen_slider_plot(gspopts);



if nargout > 0
    hout = fig;
end



function frame = get_ei_contour_frame(contours, sval)
frame = contours{sval};


function h = plot_func(frame, sval, contours, f)
xlims = getappdata(f, 'xlims');
ylims = getappdata(f, 'ylims');
zlims = getappdata(f, 'zlims');

h = plot_ei_contour(frame);
set(gca, 'XLim', xlims, 'YLim', ylims, 'ZLim', zlims);
title(sprintf('Frame %i', sval));