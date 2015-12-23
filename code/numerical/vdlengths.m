function L = vdlengths(dt, varargin)
% VDLENGTHS     Given a DelaunayTri object, calculate the length of Voronoi
%               Diagram faces between each pair of points
%
% If P is the number of points, L is a PxP matrix where each entry xi,yj,
% gives the length of the Voronoi Diagram face between point i and point j.
% If i and j do not share a face, then xi,yj is zero.
%
% Because it is convenient, this calculation use both the upper and lower
% triangles of the point-to-point arrays to do its calculation; this is
% inefficient by about 2x.
%
% Be careful that you don't have floating point imprecision or you will get
% funny triangulation around the edges.
%
% 2013-05, phli
%

opts = inputParser();
opts.addParamValue('vas', vertexAttachments(dt));
opts.addParamValue('ccs', circumcenters(dt));
opts.parse(varargin{:});
opts = opts.Results;

npoints = size(dt.X, 1);

% First get circumcenter coordinates
ccx = spalloc(npoints, npoints, dt.size(1)*3);
ccy = spalloc(npoints, npoints, dt.size(1)*3);
for i = 1:length(opts.vas)
    vas = opts.vas{i};
    triangles = dt(vas,:);
    for j = 1:length(vas)
        cc = opts.ccs(vas(j),:);
        for k = 1:3
            point = triangles(j,k);
            if point == i, continue; end
            if ccx(i,point) == 0
                ccx(i,point) = cc(1);
                ccy(i,point) = cc(2);
            elseif ccx(point,i) == 0 || ccx(point,i) == ccx(i,point);
                ccx(point,i) = cc(1);
                ccy(point,i) = cc(2);
            end
        end
    end
end

% Calculate edge lengths
ccx2 = (ccx - ccx').^2;
ccy2 = (ccy - ccy').^2;
L = sqrt(ccx2 + ccy2);