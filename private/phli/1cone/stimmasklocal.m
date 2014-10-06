function mask = stimmasklocal(conerun, rgc, stimmap, varargin)

opts = inputParser();
opts.addParamValue('az_pad_factor', 6);
opts.addParamValue('az_aspect_ratio', 1);
opts.parse(varargin{:});
opts = opts.Results;

mapsize = size(stimmap);
scale = mapsize ./ [conerun.stimulus.field_width conerun.stimulus.field_height];
bounds = autozoom_to_fit(conerun, rgc, opts.az_pad_factor, scale, opts.az_aspect_ratio);
rectmaskfunc = @(X,Y)(X > bounds(1) & X < bounds(2) & Y > bounds(3) & Y < bounds(4));
mask = funcmeshgrid(rectmaskfunc, 1:mapsize(2), 1:mapsize(1));