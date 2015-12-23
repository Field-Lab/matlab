function [in_cps_out, base_cps_out] = cpselect_pretransform(in, base, in_cps, base_cps, tform, varargin)
% CPSELECT_PRETRANSFORM    Run cpselect with a pretransformed input, keep control points in the original space
%
% usage: [in_cps, base_cps] = cpselect_pretransform(in, base, in_cps, base_cps, tform, [tform_opts])
%
% inputs: IN, BASE, IN_CPS, BASE_CPS    See documentation for CPSELECT
%         TFORM    The transform to apply to IN image before opening the
%                  cpselect interface.  IN_CPS that you pass in will also
%                  be transformed to match this space, and the IN_CPS that
%                  are returned will be properly back transformed to the
%                  original image space
%         
%         TFORM_OPTS    Remaining arguments will be passed through to
%                       IMTRANSFORM.  Particularly useful is specifying
%                       XData and YData.
%
% phli 2010-04
%


[tin, xdata, ydata] = imtransform(in, tform, varargin{:});

% IMTRANSFORM does an additional translation to put the image within bounds
% (unless you have specified XData and YData in the arguments).  So we must
% take this silent transform into account and compensate for it.
[height_in, width_in] = size(tin);
imfix_tform = cp2tform([1 1; width_in 1; 1 height_in], [xdata(1) ydata(1); xdata(2) ydata(1); xdata(1) ydata(2)], 'affine');

% Unfortunately, CPSELECT will not take empty IN_CPS/BASE_CPS so we have to
% work around with some control flow...
if isempty(in_cps) && isempty(base_cps)
    [in_cps, base_cps] = cpselect(tin, base, 'Wait', true);
elseif ~isempty(in_cps) && ~isempty(base_cps)
    % Apply both the given TFORM and the silent image fixing tform to the given IN_CPS
    in_cps = tforminv(imfix_tform, tformfwd(tform, in_cps));
    [in_cps, base_cps] = cpselect(tin, base, in_cps, base_cps, 'Wait', true);
else
    error('IN_CPS and BASE_CPS must have same number of points');
end

% Reverse the given TFORM and silent image fix tform before giving output
in_cps = tforminv(tform, tformfwd(imfix_tform, in_cps));

if nargout > 0
    in_cps_out   = in_cps;
    base_cps_out = base_cps;
end