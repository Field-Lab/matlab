function theta = angle2Points(varargin)
%ANGLE2POINTS return horizontal angle between 2 points
%
%   ALPHA = angle2Points(P1, P2),
%   Pi are either [1*2] arrays, or [N*2] arrays, in this case ALPHA is a 
%   [N*1] array. The angle computed is the horizontal angle of the line 
%   (P1 P2)
%   Result is always given in radians, between 0 and 2*pi.
%
%   See Also:
%   points2d, angle3points
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 02/03/2007.
%

%   HISTORY :

% process input arguments
if length(varargin)==3
    p1 = varargin{1};
    p2 = varargin{2};
elseif length(varargin)==1
    var = varargin{1};
    p1 = var(1,:);
    p2 = var(2,:);
end    

% ensure data have same size
if size(p1, 1)==1
    p1 = repmat(p1, [size(p2, 1) 1]);
elseif size(p2, 1)==1
    p2 = repmat(p2, [size(p1, 1) 1]);
elseif size(p1, 1)~=size(p2, 1)
    error('angle2Points: wrong size for inputs');
end

% angle line (P2 P1)
dp = p2-p1;
theta = mod(atan2(dp(:,2), dp(:,1)) + 2*pi, 2*pi);

