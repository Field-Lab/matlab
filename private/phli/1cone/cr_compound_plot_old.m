function cr_compound_plot(crstruct, rgc_indices)

crfs    = crstruct.crfs;
crrun   = crstruct.crrun;
crcells = crstruct.crcells;
crstim  = crstruct.crstim;
conerun = crstruct.conerun;
conecells = crstruct.conecells;
stabilityrun = crstruct.stablerun;

if nargin < 2 || strcmp(rgc_indices, 'all')
    rgc_indices = 1:length(crcells);
end

intensities = get_intensities(crstim);
xscales = xscalecr(crfs, crstim);

for i = 1:length(rgc_indices);
    ind = rgc_indices(i);
    cr = crfs(:,:,ind);
    if (max(cr(:)) == 0), continue; end;
    
    cones_1_and_2 = get_cones_by_weight_rank(conerun, conecells(ind), [1 2]);

    weights = conerun.cones.weights(cones_1_and_2, conecells(ind));
    weightax = subplot(length(rgc_indices), 3, sub2ind([3 length(rgc_indices)], 1, i));
    hold on;
    h1  = plot(intensities.*weights(1),   cr(:,1), 'b.');
    h2  = plot(intensities.*weights(2),   cr(:,2), 'r.');
    h12 = plot(intensities.*sum(weights), cr(:,3), 'k.');
    legend([h1 h2 h12], num2str(weights(1)), num2str(weights(2)), num2str(sum(weights)));
    drawnow;
    
    xscale = xscales{ind};
    fitax = subplot(length(rgc_indices), 3, sub2ind([3 length(rgc_indices)], 2, i));
    hold on;
    h1  = plot(intensities.*xscale(1).*sum(weights), cr(:,1), 'b.');
    h2  = plot(intensities.*xscale(2).*sum(weights), cr(:,2), 'r.');
    h12 = plot(intensities.*xscale(3).*sum(weights), cr(:,3), 'k.');
    legend([h1 h2 h12], num2str(xscale(1).*sum(weights)), num2str(xscale(2).*sum(weights)), num2str(xscale(3).*sum(weights)));
    drawnow;

    rfax = subplot(length(rgc_indices), 3, sub2ind([3 length(rgc_indices)], 3, i));
    highlight_rgba = [1 0 0; 0 1 0];
    plot_voronoi_over_rf(conerun, conecells(ind), 'highlight_cones', cones_1_and_2, 'highlight_rgba', highlight_rgba, 'mode', 'manhattan', 'title', false);
    drawnow;
end

subplot(length(rgc_indices), 3, sub2ind([3 length(rgc_indices)], 1, 1));
title('Scaled by STA cone weights')

subplot(length(rgc_indices), 3, sub2ind([3 length(rgc_indices)], 2, 1));
title('Scaled by least squares')