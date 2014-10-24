function h = drawConvexShape(pos,color,sz)

% drawConvexShape(pos,sz)
%
% draw the smooth convex hull that includes all of the points
% given by the input vector (pos, N x 2) adding space around
% the edge as determined by sz

% get the points of a circle
theta = linspace(0,2*pi,100);
x = sin(theta(:))*sz;
y = cos(theta(:))*sz;

cat_points = [];
for i=1:size(pos,1);
  tmpX = x + pos(i,1);
  tmpY = y + pos(i,2);
  cat_points = [cat_points; [tmpX,tmpY]];
end

% compute the convex hull
hull_inds = convhull(cat_points(:,1),cat_points(:,2));
hull = cat_points(hull_inds,:);

% plot the shape
hold on
h = fill(hull(:,1),hull(:,2),color);
hold off
