function varargout = points_in_stack(stack, positions, positionsy)
% POINTS_IN_STACK    Determine which of given points are within bounds of image stack
% usage: [points, points_tf, edges_xy] = points_in_stack(stack, positions_xy)
% usage: [points, points_tf, edges_x, edges_y] = points_in_stack(stack, positionsx, positionsy)
%
% inputs:   stack           Image stack struct as described in Proposal.rtf
%           positions_xy    Nx2 positions of points to check, alternatively
%                              can be given as separate vectors
%
% outputs:  points      Indices of points that are within bounds
%           points_tf   True/false vector for all points
%           edges_xy    Mx2 coordinates for line tracing the edge of
%                           bounds.  Alternatively can be returned as
%                           separate vectors.
%
% 2010-08 phli
%

if nargin < 3
    positionsx = positions(:,1);
    positionsy = positions(:,2);
else
    positionsx = positions;
end

[edgesx, edgesy] = stack_edges(stack);
[edgesxt, edgesyt] = tformfwd(stack.tforms{1}, edgesx, edgesy);
in = inpolygon(positionsx, positionsy, edgesxt, edgesyt);
inpoints = find(in);


varargout{1} = inpoints;
if nargout > 1
    varargout{2} = in;
end
if nargout == 3
    varargout{3} = [edgesxt' edgesyt'];
end
if nargout >3
    varargout{3} = edgesxt;
    varargout{4} = edgesyt;
end