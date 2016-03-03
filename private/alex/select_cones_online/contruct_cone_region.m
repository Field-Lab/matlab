function [contour, map_speck] = contruct_cone_region(radius)


x = 0;
y = 0;
th = 0:pi/250:2*pi;
xunit = round(((radius * cos(th) + x)*2))/2;
yunit = round(((radius * sin(th) + y)*2))/2;

a = [xunit; yunit]';
m = [];
for j=1:499
    if xunit(j)==xunit(j+1) && yunit(j)==yunit(j+1)
        m = [m j];
    end
end
a(m,:) = [];

[X, Y] = meshgrid(linspace(x-10,x+10, 21), linspace(y-10,y+10, 21));
in = inpolygon(X, Y, a(:,1), a(:,2));
[rows, cols] = find(in);
p = min(cols);
pp = max(cols)-p+1;

linearInd = sub2ind([pp,pp], cols-p+1, rows-p+1);
map = zeros(pp,pp);
map(linearInd) =1;

dd = imresize(map,5,'method', 'nearest');
[r, c] = find(dd,1);
contour = bwtraceboundary(dd,[r c],'W',8,Inf,'counterclockwise');
contour= round(contour/5)+0.5;
m = [];
for k = 1:length(contour)-1
    if all(contour(k,:) == contour(k+1,:))
        m = [m k];
    end
end
contour(m,:) = [];
p = round(mean(contour(:)));
contour = contour-p;

[X, Y] = meshgrid(linspace(-30,30, 61), linspace(-30,30, 61));
in = inpolygon(X, Y, contour(:,1), contour(:,2));
[r, c] = find(in);
in = in(min(r):max(r), min(c):max(c));
[r, c] = find(in);
r = r-mean(r);
c = c-mean(c);

map_speck = [r c];

