function plot_delaunay_tri(datarun, cell_spec)
% PLOT_DELAUNAY_TRI    Show the Delaunay triangulation of the receptive field centers for a given cell type
%
% usage: plot_delaunay_tri(datarun, cell_spec)
%
% 2010-01 phli
%

[cell_type_num, cell_type_name] = get_cell_type_nums(datarun, cell_spec);

delaunay = datarun.stas.delaunays{cell_type_num};
x = delaunay.coordinates(:,1);
y = delaunay.coordinates(:,2);

triplot(delaunay.triangles, x, y);
hold on;
plot(x, y, '.k');

titlestr = cell_type_name;
if isfield(delaunay, 'culled') && ~isempty(delaunay.culled)
    titlestr = [titlestr ', culled: ' delaunay.culled{1} ' ' num2str(delaunay.culled{2})];
end
title(titlestr);
        