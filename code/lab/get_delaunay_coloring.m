function datarun = get_delaunay_coloring(datarun, cell_spec)
% GET_DELAUNAY_COLORING    Datarun wrapper for DELAUNAY_COLOR
%
% usage: datarun = get_delaunay_coloring(datarun, cell_spec)
%
% 2010-02 phli
%

cell_type_num = get_cell_type_nums(datarun, cell_spec);
dt = datarun.stas.delaunays{cell_type_num};
datarun.stas.delaunays{cell_type_num}.coloring = delaunay_color(dt);