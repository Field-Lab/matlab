function tform = build_rf_tform(datarun, cell_id, tform, center_type)
% BUILD_RF_TFORM    Generate a TFORM (as from MAKETFORM) with assumed bookend translation to/from rf center
%
% usage: tform = build_rf_tform(datarun, cell_id, tform, center_type)
%
% When using affine transforms on RFs, generally makes sense to transform
% them relative to the RF center rather than relative to the image origin.
% This takes the given affine transform matrix and wraps translation
% to/from the RF center around it.
%
% 2010-02 phli
%

if nargin < 4
    center = rf_center(datarun, cell_id);
else
    center = rf_center(datarun, cell_id, center_type);
end

tform = build_affine_tform(tform);

tform = maketform('affine', [1 0 0; 0 1 0; -center(1) -center(2) 1] * tform * [1 0 0; 0 1 0; center(1) center(2) 1]);
