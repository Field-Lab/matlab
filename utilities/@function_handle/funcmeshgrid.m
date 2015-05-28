function mesh = funcmeshgrid(f, varargin)
% FUNCMESHGRID      Create meshgrid and run function on it in one step
% usage: mesh = funcmeshgrid(f, varargin)
%
% 2012-07 phli
%

% Determine whether 2D or 3D meshgrid is called for
if (length(varargin) == 3)
    gridout = cell(1,3);
else
    gridout = cell(1,2);
end

% Get mesh
[gridout{:}] = meshgrid(varargin{:});

% Run F on mesh
mesh = f(gridout{:});