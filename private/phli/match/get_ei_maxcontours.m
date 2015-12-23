function contours = get_ei_maxcontours(datarun, cell_spec, numcontours, xi, yi, interp_order, nthreads)

if nargin < 6
    interp_order = 2;
end

if nargin < 7
    nthreads = 8;
end


positions = datarun.ei.position;
[XI,YI] = meshgrid(xi,yi);

cellnums = get_cell_indices(datarun, cell_spec);
ncells = length(cellnums);
for i = 1:ncells
    cell_id = datarun.cell_ids(cellnums(i));
    if ncells > 1, fprintf('Processing cell %d (%d of %d)\n', cell_id, i, ncells); end
    
    ei = abs(get_ei_max_frame(get_ei(datarun, cell_id)));
    [cx,cy,coeffs] = compute_3dbspline_coeffs(positions(:,1), positions(:,2), ei);
    ZI = interp_3dbspline(cx,cy,coeffs,XI,YI,interp_order, nthreads);
    contours{cellnums(i)} = parse_contourc(xi, yi, ZI, numcontours);
end