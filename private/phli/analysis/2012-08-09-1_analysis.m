%% Basics
piece = '2012-08-09-1';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'udprintpath'; 'allconesprintpath'; 'keep_vars'; 'd'};


%% WN repeat rasters test
piece = '2012-08-09-1';

d.d01s = load_data(fullfile(piece, 'streamed/data001/data001'), loadopts);
d.d01s = load_sta(d.d01s, 'load_sta', []);
d.d01s = get_sta_summaries(d.d01s, {4}, 'robust_std_method', 6);
d.d01s = calc_rf_snrs(d.d01s);
[sorted, inds] = sortlownans(d.d01s.stas.medsigs);
d.d01s.stas.offM_snrsorted = d.d01s.cell_ids(inds(~isnan(sorted)));
clear sorted inds

d.d02 = load_data(fullfile(piece, 'data002'), loadopts);
cell_list_map = map_ei(d.d01s, d.d02, 'master_cell_type', {4});
d.d02.offM_highestsnr_fd01s = cell_list_map(get_cell_indices(d.d01s, d.d01s.stas.offM_snrsorted(end:-1:end-10)));
clear cell_list_map;

% Check trigger spacing
plot(diff(d.d02.triggers), 'o-');

for i = 1:length(d.d02.offM_highestsnr_fd01s)
    cellid = d.d02.offM_highestsnr_fd01s{i};
    if isempty(cellid), continue; end
    
    figure();
    rasterphli(d.d02, cellid, d.d02.triggers(1:6:end));
    axis ij;
end
clear i cellid