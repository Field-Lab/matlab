function varargout = drawPoint(varargin)
%DRAWPOINT draw the point on the axis.
%
%   drawPoint(X, Y);
%   Draws points defined by coordinates X and Y.
%   X and Y are N*1 array, with N being number of points to be drawn.
%   If coordinates of points lie outside the visible area, points are
%   not drawn.
%
%   drawPoint(COORD);
%   Packs coordinates in a single [N*2] array.
%
%   drawPoint(..., OPT);
%   Draws each point with given option. OPT is a series of arguments pairs
%   compatible with 'plot' model.
%
%
%   H = drawPoint(...) also return a handle to each of the drawn points.
%
%   See also
%   points2d, clipPoints
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY
%   23/02/2004: add more documentation. Manage different kind of inputs. 
%     Does not draw points outside visible area.
%   26/02/2007: update processing of input arguments.


% process input arguments
var = varargin{1};
if size(var, 2)==1
    % points stored in separate arrays
    px = varargin{1};
    py = varargin{2};
    varargin(1:2) = [];
else
    % points packed in one array
    px = var(:, 1);
    py = var(:, 2);
    varargin(1) = [];
end

% ensure we have column vectors
px = px(:);
py = py(:);

% default drawing options, but keep specified options if it has the form of
% a bundled string
if length(varargin)~=1
    varargin = {'linestyle', 'none', 'marker', 'o', 'color', 'b', varargin{:}};
end

% limits of the current axis
lim = get(gca, 'xlim');
xmin = lim(1); 
xmax = lim(2);
lim = get(gca, 'ylim');
ymin = lim(1);
ymax = lim(2);

% compute indices of points inside visible area
ok = px>=xmin;
ok = ok & px<=xmax;
ok = ok & py>=ymin;
ok = ok & py<=ymax;

% plot the points, using specified drawing options
h = plot(px(ok), py(ok), varargin{:});

% process output arguments
if nargout>0
    varargout{1}=h;
end
