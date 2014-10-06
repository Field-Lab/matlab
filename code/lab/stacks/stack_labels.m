function [cells, rows] = stack_labels()
% STACK_LABELS  Data for translating named indices into numeric
%
% 2010-09 phli
%

cells = struct();
cells.array       = {1,1};
cells.array_image = {1,2};
cells.alive_montage_lores = {2,1};
cells.alive_montage_hires = {2,2};
cells.photographic_mapping_array_edges = {4,1};
cells.photographic_mapping_array       = {4,2};
cells.photographic_mapping_pm32        = {4,3};
cells.photographic_mapping_pm10        = {4,4};
cells.photographic_mapping_pm2         = {4,5};

if nargout > 1
    rows = struct();
    rows.array                = 1;
    rows.alive_full           = 2;
    rows.alive_montage        = 2;
    rows.alive_partial        = 3;
    rows.photographic_mapping = 4;
    rows.fixed_lores          = 10;
    rows.fixed_hires          = 11;
end