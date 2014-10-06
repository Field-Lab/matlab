function points = traced_cell_points(soma,axon)
% traced_cell_points     take a soma location and radius, and axon trajectory, and
%                           create a list of points which can be plotted as a line
%
% usage:  points = traced_cell_points(soma,axon)
%
% arguments:     soma - 2x2 matrix, first row is coordinates of soma center, second row is where axon leaves soma
%                           assumes a circular soma
%                axon - Ax2 matrix, coordinates of the axon trajectory, starting at the edge of the soma
%
% outputs:     points - Nx2 matrix of coordinates
%
%
% example: points = traced_cell_points(axons{axon_id}(1:2,:),axons{axon_id}(2:end,:));
%
% 2010-03  gauthier
%

% ToDo-PHL-2010-05-24: Allow to take single argument axon trace, assuming
% first point is middle of cell body, second point is where exits cell
% body.


% get angle, diameter in polar coords
[theta, rho] = cart2pol(soma(2,1)-soma(1,1),soma(2,2)-soma(1,2));

% make offsets
offsets = [0:0.1:2*pi 2*pi];

% compute soma points in polar coords
[soma_x,soma_y] = pol2cart(theta+offsets', repmat(rho,length(offsets),1));

% concatenate
points = [soma_x+soma(1,1) soma_y+soma(1,2); axon];

