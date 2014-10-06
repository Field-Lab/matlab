function datarun = simplify_rf_contours(datarun, cell_spec, contour_indices)
% SIMPLIFY_RF_CONTOURS    Reduce each rf's contour polygons to the single largest area polygon
% usage: datarun = simplify_rf_contours(datarun, cell_spec, contour_indices)
%
% Saves results to DATARUN.STAS.RF_CONTOURS_SIMPLE
%
% phli 2010-03
%


cell_nums = get_cell_indices(datarun, cell_spec);
for i = 1:length(cell_nums)
    cell_num = cell_nums(i);

    if nargin < 3
        contour_indices = 1:length(datarun.stas.rf_contours{cell_num});
    end

    for contour_num = contour_indices
        contour = datarun.stas.rf_contours{cell_num}{contour_num};
        datarun.stas.rf_contours_simple{cell_num}{contour_num} = simplify_contour(contour);
    end
end