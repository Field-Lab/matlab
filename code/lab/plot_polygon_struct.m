function plot_polygon_struct(polygon_struct, varargin)
% PLOT_POLYGON_STRUCT    Display a polygon struct
%
% usage: plot_polygon_struct(polygon_struct, facecolor, bgcolor, alpha)
%
% Polygon struct is the format used by PolygonClip and many SNL-E analyses.
% Basically has fields 'x', 'y', 'hole'.  
%
% I have not discovered a good way to handle plotting polygons with holes
% in Matlab.  Here I took the simplistic approach of plotting holes as
% secondary polygons filled in with the background color, but this may be
% wrong for some complicated situations with multiple polygon_structs
% overlapping each other.
%
% 2010-02 phli
%

opts = inputParser;
opts.addParamValue('facecolor', 'b');
opts.addParamValue('bgcolor', 'w');
opts.addParamValue('alpha', 1);
opts.addParamValue('linecolor', 'k');
opts.addParamValue('linewidth', 1);
opts.parse(varargin{:});
opts = opts.Results;

hold on;
for i = 1:numel(polygon_struct);
    if polygon_struct(i).hole == 0
        % Have to hack the face color a little to allow 'none'
        fill(polygon_struct(i).x, polygon_struct(i).y, 'k', 'FaceAlpha', opts.alpha, 'EdgeColor', opts.linecolor, 'FaceColor', opts.facecolor, 'LineWidth', opts.linewidth);
    end
end

for i = 1:numel(polygon_struct);
    if polygon_struct(i).hole == 1
        fill(polygon_struct(i).x, polygon_struct(i).y, opts.bgcolor, 'EdgeColor', opts.linecolor, 'LineWidth', opts.linewidth);
    end
end