function contours = compute_ei_contours(ei, positions, numcontours, xi, yi, interp_order, nthreads)

if nargin < 6
    interp_order = 2;
end

if nargin < 7
    nthreads = 8;
end


[XI,YI] = meshgrid(xi,yi);


eiframes = size(ei,2);
contours = cell(eiframes,1);
for i = 1:eiframes
    if eiframes > 1, fprintf('Processing frame %d of %d\n', i, eiframes); end
    
    [cx,cy,coeffs] = compute_3dbspline_coeffs(positions(:,1), positions(:,2), ei(:,i));
    ZI = interp_3dbspline(cx,cy,coeffs,XI,YI,interp_order, nthreads);
    contours{i} = parse_contourc(xi, yi, ZI, numcontours);
end