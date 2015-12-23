function plot_parasol_contours(datarun)

hold on
ids = datarun.cell_types{2}.cell_ids;

ai = parse_stack_index('array');
array_tform_inv = datarun.stacks{11,3}.tforms_inv{ai{:}};

for i = 1:length(ids)
    ei = get_ei(datarun, ids(i));
    max_frame = get_ei_max_frame(ei);
    contours = hex_contour(datarun.ei.position(:,1), datarun.ei.position(:,2), max_frame, 10);
    plot_contours(contours, array_tform_inv);
end



function plot_contours(contours, tform)

for i = 1
    contour = contours(i);
    
    for j = 1:length(contour.paths)
        path = contour.paths{j};
        if nargin > 1, path = tformfwd(tform, path); end
        h = plot(path(:,1),path(:,2), 'b--');
        set(h, 'Color', [0.5 0.5 1]);
        set(h, 'LineWidth', 1);
    end
end