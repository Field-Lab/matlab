function [XI, YI, ZI] = interp_scattered_stack_points(stack, points, varargin)

opts = inputParser;
opts.addParamValue('xlim', 1:size(stack, 1));
opts.addParamValue('ylim', 1:size(stack, 2));
opts.addParamValue('method', 'nearest');
opts.addParamValue('verbose', true);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
opts = opts.Results;



[X,V] = stack_points_to_XV(points);

if opts.verbose; disp('Building scattered point interpolator'); end
TSI = TriScatteredInterp(X, V, opts.method);

if opts.verbose; disp('Building X-Y meshgrid'); end
[XI, YI] = meshgrid(opts.ylim, opts.xlim);



if opts.verbose; fprintf('Interpolating Z values from scattered points; this can take a while...'); end
ZI = TSI(XI, YI);
if opts.verbose; disp(' Done'); end



if ~strcmp(opts.method, 'nearest') && any(isnan(ZI(:)))
    if opts.verbose; fprintf('Using nearest neighbor to fill in locations with no linear interpolant; this can take a while...'); end
    TSI2 = TriScatteredInterp(X, V, 'nearest');
    ZI2 = TSI2(XI, YI);
    ZI(isnan(ZI)) = ZI2(isnan(ZI));
    if opts.verbose; disp(' Done'); end
end