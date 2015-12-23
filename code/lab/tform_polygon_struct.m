function pgs = tform_polygon_struct(pgs, tform)
% TFORM_POLYGON_STRUCT    Apply a Matlab 2D transformation to a polygon struct object
% usage: pgs = tform_polygon_struct(pgs, tform)
%
% phli 2010-04
%

for i = 1:numel(pgs)
    u = pgs(i).x;
	v = pgs(i).y;

	[x, y] = tformfwd(tform, u, v);
	
	pgs(i).x = x;
	pgs(i).y = y;
end
