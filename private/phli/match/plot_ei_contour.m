function hout = plot_ei_contour(contour)

oldhold = ishold;
if ~oldhold
    cla
end
hold('on');

h = gcf;


for i = 1:length(contour)
    path = contour(i).path;
    elevation = contour(i).elevation;
    if elevation <= 0
        color = 'b';
    else
        color = 'r';
    end
    plot3(path(1,:), path(2,:), repmat(elevation, size(path(1,:))), color);
end


if ~oldhold
    hold off
end

if nargout > 0
    hout = h;
end