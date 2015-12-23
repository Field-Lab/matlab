function new_points = project_to_orthogonal_susbpace(points, vector)
% project_to_orthogonal_susbpace     project points to a subspace orthogonal to a given vector
%
% usage:  new_points = project_to_orthogonal_susbpace(points, vector)
%
% arguments:     points - MxN matrix, each row a point in N-d space
%                vector - N vector
%
% outputs:     new_points - MxN matrix, each row a point in N-d space
%
%
% 2009-10 gauthier
%


% ensure column
vector = reshape(vector,[],1);

% normalize
vector = vector / norm(vector);

% compute projection
new_points = points - (points * vector) * vector';


% plot first 3 dimensions
if 0

    figure(1);clf;hold on

    % plot old points, connecting to new points
    for pp=1:size(points,1)
        plot3([points(pp,1) new_points(pp,1)],[points(pp,2) new_points(pp,2)],[points(pp,3) new_points(pp,3)],'b')
    end

    % plot new points
    plot3(new_points(:,1),new_points(:,2),new_points(:,3),'.k')

    % plot vector
    plot3([0 vector(1)],[0 vector(2)],[0 vector(3)],'r')

    axis equal

end
