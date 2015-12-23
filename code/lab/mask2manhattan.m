function [edges_x, edges_y] = mask2manhattan(mask)
% usage: [edges_x, edges_y] = mask2manhattan(mask)
% 
% Example: mask = poly2mask(x, y, m, n);
%          [edges_x, edges_y] = mask2manhattan(mask);
%          plot(edges_x, edges_y);
%
% 2011-07 phli

% Find on pixels, row and column indices
[is,js] = find(mask);

edges_x = cell(1, length(is));
edges_y = cell(1, length(js));
for k = 1:length(is)
    i = is(k);
    j = js(k);
    [edges_x{k}, edges_y{k}] = get_unneighbored_edges(mask, i, j);
end
edges_x = cell2mat(edges_x);
edges_y = cell2mat(edges_y);


% For the given pixel (i,j), which should be an on pixel, check whether the
% neighboring pixels in all 4 directions are also on.  If they are on, then
% we don't need an edge on that side, but if they aren't, then add that
% edge.
function [xes, ys] = get_unneighbored_edges(mask, i, j)
[h,w] = size(mask);

left   = j-0.5;
right  = j+0.5;
top    = i-0.5;
bottom = i+0.5;

xes = [];
ys  = [];
if j > 1 && ~mask(i,j-1)
    leftedgex = [left; left];
    leftedgey = [top;  bottom];
    xes(1:2, end+1) = leftedgex;
    ys( 1:2, end+1) = leftedgey;
end
if j < w && ~mask(i,j+1)
    rightedgex = [right; right];
    rightedgey = [top;   bottom];
    xes(1:2, end+1) = rightedgex;
    ys( 1:2, end+1)  = rightedgey;
end
if i > 1 && ~mask(i-1,j)
    topedgex = [left; right];
    topedgey = [top;  top];
    xes(1:2, end+1) = topedgex;
    ys( 1:2, end+1) = topedgey;    
end
if i < h && ~mask(i+1,j)
    bottomedgex = [left;   right];
    bottomedgey = [bottom; bottom];
    xes(1:2, end+1) = bottomedgex;
    ys( 1:2, end+1) = bottomedgey;
end