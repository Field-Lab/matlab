function ax = plot_rf_stimmap(datarun, cellid, map, varargin)
% PLOT_RF_STIMMAP   Plot STA RF with stimulus map regions overlaid
% usage: ax = plot_rf_stimmap(datarun, cellid, map, varargin)
%
% Can accept maps as either Manhattan line segments (e.g. output from
% MAP2MANHATTAN), or as polylines (e.g. output of
% SEGS2POLY(MAP2MANHATTAN(...))).
%
% 2012-10 phli
%

opts = inputParser();
opts.KeepUnmatched = true;
opts.addParamValue('bounds', []);
opts.addParamValue('autozoom', true);
opts.addParamValue('az_pad_factor', 4);
opts.addParamValue('az_aspect_ratio', 1);
opts.addParamValue('marksthresh', 3.5);
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;

if ~isempty(opts.bounds), opts.autozoom = false; end

rfopts = keepfields(unmatched, 'foa', 'coordinates', 'com', 'sig_stix', 'polarity', 'color_transform', 'array', 'fit', 'title', 'pdtitle', 'ticks', '-force');
if isfield(unmatched, 'scaled_up'), rfopts.scale = unmatched.scaled_up; end
datarun = setsimplemarks(datarun, cellid, opts.marksthresh);
ax = plot_rf(datarun, cellid, rfopts);
hold on;


scalemap = [1 1] ./ [datarun.stimulus.stixel_width datarun.stimulus.stixel_height];
patchopts = keepfields(unmatched, 'colors', 'cmf', 'fill', 'fillcolors', 'patchopts', '-force');
patchopts.scale = scalemap;

% Can accept either line segments or polyline (see help)
if isstruct(map{1})
    patchmanhattan(map, patchopts);
else
    patchpolylines(map, patchopts);
end


if opts.autozoom
    opts.bounds = autozoom_to_fit(datarun, cellid, opts.az_pad_factor, 1, opts.az_aspect_ratio);
end
if ~isempty(opts.bounds)
    axis(opts.bounds);
end