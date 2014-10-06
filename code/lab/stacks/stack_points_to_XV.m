function [X,V] = stack_points_to_XV(points)
% Convert XYZ points from stack_point_picker GUI format to format expected
% by TriScatteredInterp.
X = vertcat(points{:});
V = [];
for i = 1:length(points)
    section_points = points{i};
    for j = 1:size(section_points, 1)
        V(end+1,1) = i;
    end
end