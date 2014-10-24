function M = getMixMat(locs,params)

d = getDistances(locs,locs);
c = size(d,1);

switch params.type
    case 'dist'
        M = params.p(1)*exp(-(d/params.p(2)).^params.p(3));
    case 'e-neighborhood'
        M = double(d<params.p(1))*params.p(2);
    case 'k-nearest'
        k = kNearestNeighborGraph(max(d(:))-d,params.p(1));
        M = double(k>0)*params.p(2);
end

M = M - diag(diag(M));
M = M + diag(ones(1,c));
        