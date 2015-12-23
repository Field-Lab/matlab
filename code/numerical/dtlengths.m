function L = dtlengths(dt)
% DTLENGTHS     Given a DelaunayTri object, calculate the length of edges
%               between each pair of points
%
% If P is the number of points, L is a PxP matrix where each entry xi,yj,
% gives the length of the edge between point i and point j.  If i and j do 
% not share an edge, then xi,yj is Inf.
%
% Be careful that you don't have floating point imprecision or you will get
% funny triangulation around the edges.
%
% 2013-05, phli
%

npoints = size(dt.X, 1);
L = inf(npoints, npoints);
for i = 1:size(dt.Triangulation, 1)
    tri = dt.Triangulation(i,:);
    pairs = tri([1 2; 2 3; 3 1]);
    
    for j = 1:3
        p1 = pairs(j,1); 
        p2 = pairs(j,2);
        if L(p1,p2) == Inf
            L(p1,p2) = sqrt(sum(diff(dt.X([p1 p2],:)).^2));
            L(p2,p1) = L(p1,p2);
        end
    end
end