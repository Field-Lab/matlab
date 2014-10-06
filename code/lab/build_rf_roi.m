function datarun = build_rf_roi(datarun, cell_spec, do_plot)
% BUILD_RF_ROI    Create a Region of Interest for RF calculations
%
% usage: datarun = build_rf_roi(datarun, cell_spec, [do_plot])
%
% For mosaic calculations, it is useful to have an roi indicating what area
% the mosaic covers.
%

if nargin < 3
    do_plot = false;
end

cell_type_num = get_cell_type_nums(datarun, cell_spec);
delaunay = datarun.stas.delaunays{cell_type_num};
polygons = triangles_to_polygons(delaunay.triangles, delaunay.coordinates);


% Initialize the ROI as the first triangle
roi = polygons(1);

% Union the ROI with every triangle
union_mode = 3; % Arcane setting for PolygonClip
for i = 2:numel(polygons)
    roi = PolygonClip(roi, polygons(i), union_mode);
end

datarun.stas.delaunays{cell_type_num}.roi = roi;

if do_plot
    figure;
    plot_polygon_struct(roi);
end