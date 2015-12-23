function resliced = reslice_mono_stack(stack, XI, YI, ZI, varargin)

opts = inputParser;
opts.addParamValue('method', []);
opts.addParamValue('verbose', true);
%opts.addParamValue('thickness', 1);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
opts = opts.Results;

% % For thickness opt, but not yet successfully implemented
%range = 1:opts.thickness;
%mid = (opts.thickness - 1) / 2 + 1;
%offsets = range - round(mid);

switch opts.method
    case 'nearest'
        if opts.verbose; disp('Picking points from interpolated Z values'); end
        I = sub2ind(size(stack), YI, XI, round(ZI)); % Not sure why this is Y first; has to do with meshgrid versus ndgrid, but trying to switch to ndgrid gave other weirdness
        resliced = stack(round(I));

    case 'linear'
        stack = double(stack);
        if opts.verbose; fprintf('Linearly interpolating image from Z values; this can take a while...'); end
        resliced = interp3(stack, XI, YI, ZI, 'linear');
        if opts.verbose; disp(' Done'); end
end