function joined = segs2poly(segs)
% SEGS2POLY     Recursively join line segments into polygons
% usage: joined = segs2poly(segs)
%
% Format for segs should be struct with x and y fields, as output by
% MAP2MANHATTAN.  MAP2MANHATTAN will give a cell array of these structs,
% one for each map index, so if you want to send output from MAP2MANHATTAN
% directly to here, use the wrapper NYC2POLY.
%
% This will only work for cleanly separated regions without holes.
% Otherwise you will drive Computer crazy!!
%
% Output from this can be sent to PATCHPOLYLINES.
%
% See also: MAP2MANHATTAN, NYC2POLY, PATCHPOLYLINES
%
% 2012-06 phli
%


hash = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:length(segs.x)
    for p = 1:2
        x = segs.x(p,i);
        y = segs.y(p,i);
        
        if ~hash.isKey(x)
            hash(x) = containers.Map('KeyType', 'double', 'ValueType', 'any');
        end
        
        xhash = hash(x);
        if ~xhash.isKey(y)
            xhash(y) = [];
        end
        xhash(y) = [xhash(y) i];
    end
end

visited = false(length(segs.x), 1);
joined = {};
while(any(~visited))
    tovisit = find(~visited);
    tovisit = tovisit(1);
    visited(tovisit) = true;
    [prop visited] = propagate(segs, hash, visited, segs.x(2,tovisit), segs.y(2,tovisit));
    [back visited] = backprop( segs, hash, visited, segs.x(1,tovisit), segs.y(1,tovisit));
    joined{end+1} = [back, prop];
end


% Recursive
function [prop visited] = propagate(segs, hash, visited, startx, starty)
prop = [startx; starty];
next = get_next(startx, starty, hash, visited);

% End condition
if isempty(next), return; end

% Recursion
next = next(1);
visited(next) = true;
[nextx nexty] = get_nextxy(startx, starty, next, segs);
[nextprop visited] = propagate(segs, hash, visited, nextx, nexty);
prop = [prop nextprop];


function [back visited] = backprop(segs, hash, visited, startx, starty)
back = [startx; starty];
next = get_next(startx, starty, hash, visited);

% End condition
if isempty(next), return; end

% Recursion
next = next(1);
visited(next) = true;
[nextx nexty] = get_nextxy(startx, starty, next, segs);
[nextback visited] = backprop(segs, hash, visited, nextx, nexty);
back = [nextback back];


function connected = get_connected(x,y,hash)
hashx = hash(x);
connected = hashx(y);

function unvisited = get_unvisited(tocheck, visited)
unvisited = tocheck(~visited(tocheck));

function next = get_next(x, y, hash, visited)
connected = get_connected(x, y, hash);
next = get_unvisited(connected, visited);

function [nextx nexty] = get_nextxy(prevx, prevy, nextseg, segs)
xys = [segs.x(:,nextseg) segs.y(:,nextseg)];
use = find(~all(xys' == [prevx prevy; prevx prevy]'));
nextx = xys(use,1);
nexty = xys(use,2);
