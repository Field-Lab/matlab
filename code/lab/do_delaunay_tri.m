function datarun = do_delaunay_tri(datarun, cell_spec, do_plot)
% DO_DELAUNAY_TRI Do the Delaunay triangulation of receptive field centers for a given cell type
%
% usage: datarun = do_delaunay_tri(datarun, cell_spec, [do_plot])
%
% 2010-01 phli
%

if nargin < 3
    do_plot = false;
end

[cell_nums, cell_type, cell_type_num] = get_cell_indices(datarun, cell_spec);  %#ok<ASGLU> cell_type is not used
coords = vertcat(datarun.stas.rf_coms{cell_nums});
datarun.stas.delaunays{cell_type_num} = struct('triangles', delaunay(coords(:,1), coords(:,2)),'coordinates', coords);

if nargout < 1 || do_plot
    plot_delaunay_tri(datarun, cell_type_num)
end