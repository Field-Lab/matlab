function points = svg_path_to_coordinates(path_text)
% svg_path_to_coordinates     Convert a string of path coordinates in SVG file to matrix
%
% usage:  points = svg_path_to_coordinates(path_text)
%
% arguments:     path_text - string of path coordinates from SVG file
%
% outputs:     points - Nx2 matrix
%
% http://www.w3.org/TR/SVG11/paths.html#PathData
%
%
% 2009-11  gauthier
%


% parse out the text for each point
point_texts = textscan(path_text,'%s');
point_texts = point_texts{1};

% initialize
points = [];

% for each one...
for pp=1:length(point_texts)
    point_text = point_texts{pp};
    
    % identify what kind of point it is
    switch point_text(1)
        
        case {'M','L'} % ordinary point
            
            % get the x- and y-coordinate, and save to points
            temp = textscan(point_text(2:end),'%f','delimiter',',');
            points(size(points,1)+1,:) = temp{1}';
            
       case 'C' % bezier curve
            
            % save only the last two numbers as the x- and y-coordinate
            temp = textscan(point_text(2:end),'%f','delimiter',',');
            points(size(points,1)+1,:) = temp{1}(end-1:end)';
            
        otherwise
            error('Unexpected type of point ''%s''!',point_text(1))
            
    end
end

