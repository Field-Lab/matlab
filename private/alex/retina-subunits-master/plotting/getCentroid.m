function centroid = getCentroid(locs)

tightBoundMn = min(locs);
tightBoundMx = max(locs);
boundX = [tightBoundMn(1) tightBoundMx(1)];
boundY = [tightBoundMn(2) tightBoundMx(2)];

rng = mean(boundX);
boundX = (boundX - rng)*0.2;
boundX = boundX + rng;

rng = mean(boundY);
boundY = (boundY - rng)*0.2;
boundY = boundY + rng;

xInd = linspace(boundX(1),boundX(1),20);
yInd = linspace(boundY(2),boundY(2),20);
[X Y] = meshgrid(xInd,yInd);
X = X(:);
Y = Y(:);
foo = getDistances(locs,[X Y]);
[~, ind] = max(mean(foo));
centroid = [X(ind) Y(ind)];

centroid = mean(locs,1);