function masks = make_voronoi_masks(V, C, width, height, scale)
% MAKE_VORONOI_MASKS    Convert Voronoi cone tesselation to logical pixel masks
% usage: masks = make_voronoi_masks(V, C, width, height)
%
% inputs: V, C              from VORONOIN
%         WIDTH, HEIGHT     mask width and height
%
% outputs: MASKS    Cell array of length C, each element a height-x-width,
%                   sparse logical array
%
% Takes polygon data from Voronoi tesselation and converts to a stack of
% logical pixel masks.
%
% See also: VORONOIN
%
% phli 2011-01
%

if nargin < 5
    scale = 1;
end


masks = cell(length(C), 1);
for i = 1:length(C)
    % First get all the corners
    corners = C{i};    
    xpoints = V(corners, 1);
    ypoints = V(corners, 2);

    % Filter out corners that are out of bounds
    inbounds =            xpoints > 0    & ypoints > 0;
    inbounds = inbounds & xpoints ~= Inf & ypoints ~= Inf;
    corners = corners(inbounds);
 
    % Do we have enough corners to make a polygon?
    if length(corners) < 3, continue; end
    
    xpoints = V(corners, 1);
    ypoints = V(corners, 2);
    
    % Need to do 0.5 shifts before/after rescaling to account for pixel boundaries starting at 0.5 in Matlab's scheme
    % Corrected 2012-02, so scaled up masks from before this date must be off by 1/2 pixel!
    xpoints = (xpoints - 0.5) .* scale + 0.5;
    ypoints = (ypoints - 0.5) .* scale + 0.5;
    
    masks{i} = sparse(poly2mask(xpoints, ypoints, height, width));
end