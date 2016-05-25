function [RGB_weights, movie] = RGB_to_BW(datarun, cell, varargin)
% cell is the cell index, NOT the vision id
% if you don't have the index, run
% cell  = find(datarun.cell_ids == cell_id); first
% if a color movie is input, the function also returns the psuedoBW movie

p = inputParser;
p.addParameter('color_movie', 0)
p.parse(varargin{:});
color_movie = p.Results.color_movie;
clear p

% Error if the right part of datarun doesn't exist
try 
    isempty(datarun.vision.timecourses(cell).r);
catch
    error('Requires datarun.vision.timecourses')
end

% make a matrix of the time courses
RGB_matrix=[datarun.vision.timecourses(cell).r ...
    datarun.vision.timecourses(cell).g ...
    datarun.vision.timecourses(cell).b];

% SVD
[U,S,V]=svd(RGB_matrix);

% Take the first component
RGB_weights=V(:,1);

% Warning if the first component doesn't dominate
if S(1,1)/S(2,2) < 10
    warning(['Ratio of first principle component to second is only ' num2str(S(1,1)/S(2,2)) '. The psuedoBW movie might not be accurate.'])
end

if any(color_movie ~= 0)
    % Turn RGB movie into greyscale movie
    movie=squeeze(RGB_weights(1)*color_movie(:,:,1,:)+ ...
        RGB_weights(2)*color_movie(:,:,2,:)+ ...
        RGB_weights(3)*color_movie(:,:,3,:));
end

end
