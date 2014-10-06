function parsed = parse_contourc(varargin)

if nargin == 1 && size(varargin{1}, 1) == 2
    % Assume that contour array has already been computed
    C = varargin{1};
else
    C = contourc(varargin{:});
end

parsed = struct([]);
colnum = 1;
contournum = 1;
numcols = size(C, 2);
while colnum < numcols
    parsed(contournum).elevation = C(1,colnum);
    npoints = C(2,colnum);
    parsed(contournum).path = C(:,(colnum+1):(colnum+npoints));
    
    colnum = colnum + npoints + 1;
    contournum = contournum + 1;
end