%%

piece = '2012-08-21-0';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d1 = load_data(fullfile('/snle/acquisition/', piece, '/data001/data001'), loadopts);
d2 = load_data(fullfile('/snle/acquisition/', piece, '/data002/data002'), loadopts);

d1 = load_sta(d1, 'load_sta', []);
d1 = get_sta_summaries(d1, {4}, 'robust_std_method', 6);
d1 = calc_rf_snrs(d1);
[sorted, inds] = sortlownans(d1.stas.maxsigs);
d1.stas.offM_snrsorted = d1.cell_ids(inds(~isnan(sorted)));
clear sorted inds

cell_list_map = map_ei(d1, d2, 'master_cell_type', {4});
d2.offM_highestsnr_fd01s = cell_list_map(get_cell_indices(d1, d1.stas.offM_snrsorted(end:-1:end-10)));
clear cell_list_map;

% Check trigger spacing
plot(diff(d2.triggers), 'o-');

for i = 1:length(d2.offM_highestsnr_fd01s)
    cellid = d2.offM_highestsnr_fd01s{i};
    if isempty(cellid), continue; end
    
    figure();
    rasterphli(d2, cellid, d2.triggers(1:6:end));
    axis ij;
end
clear i cellid


%%
piece = '2012-08-21-2';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d1 = load_data(fullfile('/snle/acquisition/', piece, '/data000/data000'), loadopts);
d2 = load_data(fullfile('/snle/acquisition/', piece, '/data002/data002'), loadopts);

d1 = load_sta(d1, 'load_sta', []);
d1 = get_sta_summaries(d1, {4}, 'robust_std_method', 6);
d1 = calc_rf_snrs(d1);
[sorted, inds] = sortlownans(d1.stas.maxsigs);
d1.stas.offM_snrsorted = d1.cell_ids(inds(~isnan(sorted)));
clear sorted inds

cell_list_map = map_ei(d1, d2, 'master_cell_type', {4});
d2.offM_highestsnr_fd01s = cell_list_map(get_cell_indices(d1, d1.stas.offM_snrsorted(end:-1:end-20)));
clear cell_list_map;

% Check trigger spacing
plot(diff(d2.triggers), 'o-');

for i = 1:length(d2.offM_highestsnr_fd01s)
    cellid = d2.offM_highestsnr_fd01s{i};
    if isempty(cellid), continue; end
    
    figure();
    rasterphli(d2, cellid, d2.triggers(1:6:end));
    axis ij;
end
clear i cellid


%%
piece = '2012-08-21-2';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d1 = load_data(fullfile('/snle/acquisition/', piece, '/data001/data001'), loadopts);
d2 = load_data(fullfile('/snle/acquisition/', piece, '/data003/data003'), loadopts);

d1.offM = [153
    196
    391
    408
    887
    1268
    1777
    2223
    2266
    2708
    2828
    3033
    4956
    5447
    6005
    6721
    7682]';

cell_list_map = map_ei(d1, d2, 'master_cell_type', {4});
d2.offM_f01 = cell_list_map(get_cell_indices(d1, d1.offM));
clear cell_list_map;

% Check trigger spacing
plot(diff(d2.triggers), 'o-');

for i = 1:length(d2.offM_f01)
    cellid = d2.offM_f01{i};
    if isempty(cellid), continue; end
    
    figure();
    rasterphli(d2, cellid, d2.triggers(1:2:end));
    axis ij;
end
clear i cellid

%%
piece = '2012-09-13-2';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d1 = load_data(fullfile('/snle/acquisition/', piece, '/data001/data001'), loadopts);
d2 = load_data(fullfile('/snle/acquisition/', piece, '/data002/data002'), loadopts);

d1.onM = [1909 2225 2267 2449 2795 4056 4337 5556 6079 6153 6409 6679 6802 6917 7292];
d1.rgcs = d1.onM;

cell_list_map = map_ei(d1, d2);
d2.rgcs = cell_list_map(get_cell_indices(d1, d1.rgcs));
d2.rgcs
clear cell_list_map;

% Check trigger spacing
plot(diff(d2.triggers), 'o-');

for i = 1:length(d2.rgcs)
    cellid = d2.rgcs{i};
    if isempty(cellid), continue; end
    
    figure();
    rasterphli(d2, cellid, d2.triggers(1:6:end));
    axis ij;
end
clear i cellid


%%
piece = '2012-09-18-0';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d1 = load_data(fullfile('/snle/acquisition/', piece, '/data001/data001'), loadopts);
d2 = load_data(fullfile('/snle/analysis/',    piece, '/data002/data002'), loadopts);

d1 = load_sta(d1, 'load_sta', []);
d1 = get_sta_summaries(d1, {4}, 'robust_std_method', 6);
d1 = calc_rf_snrs(d1);
[sorted, inds] = sortlownans(d1.stas.maxsigs);
d1.stas.offM_snrsorted = d1.cell_ids(inds(~isnan(sorted)));
clear sorted inds

cell_list_map = map_ei(d1, d2, 'master_cell_type', {4});
d2.offM_highestsnr_fd01s = cell_list_map(get_cell_indices(d1, d1.stas.offM_snrsorted(end:-1:end-20)));
clear cell_list_map;

% Check trigger spacing
plot(diff(d2.triggers), 'o-');

for i = 1:length(d2.offM_highestsnr_fd01s)
    cellid = d2.offM_highestsnr_fd01s{i};
    if isempty(cellid), continue; end
    
    figure();
    rasterphli(d2, cellid, d2.triggers(1:6:end));
    axis ij;
end
clear i cellid


%%
piece = '2012-09-18-3';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d1 = load_data(fullfile('/snle/acquisition/', piece, '/data001/data001'), loadopts);
d2 = load_data(fullfile('/snle/acquisition/', piece, '/data002/data002'), loadopts);

d1 = load_sta(d1, 'load_sta', []);
d1 = get_sta_summaries(d1, {4}, 'robust_std_method', 6);
d1 = calc_rf_snrs(d1);
[sorted, inds] = sortlownans(d1.stas.maxsigs);
d1.stas.offM_snrsorted = d1.cell_ids(inds(~isnan(sorted)));
clear sorted inds

cell_list_map = map_ei(d1, d2, 'master_cell_type', {4});
d2.offM_highestsnr_fd01s = cell_list_map(get_cell_indices(d1, d1.stas.offM_snrsorted(end:-1:end-20)));
clear cell_list_map;

% Check trigger spacing
plot(diff(d2.triggers), 'o-');

for i = 1:length(d2.offM_highestsnr_fd01s)
    cellid = d2.offM_highestsnr_fd01s{i};
    if isempty(cellid), continue; end
    
    figure();
    rasterphli(d2, cellid, d2.triggers(1:6:end));
    axis ij;
end
clear i cellid


%%
piece = '2012-09-18-3';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d1 = load_data(fullfile('/snle/acquisition/', piece, '/data003/data003'), loadopts);
d2 = load_data(fullfile('/snle/acquisition/', piece, '/data004/data004'), loadopts);

d1 = load_sta(d1, 'load_sta', []);
d1 = get_sta_summaries(d1, {4}, 'robust_std_method', 6);
d1 = calc_rf_snrs(d1);
[sorted, inds] = sortlownans(d1.stas.maxsigs);
d1.stas.offM_snrsorted = d1.cell_ids(inds(~isnan(sorted)));
clear sorted inds

cell_list_map = map_ei(d1, d2, 'master_cell_type', {4});
d2.offM_highestsnr_fd01s = cell_list_map(get_cell_indices(d1, d1.stas.offM_snrsorted(end:-1:end-20)));
clear cell_list_map;

% Check trigger spacing
plot(diff(d2.triggers), 'o-');

f = figure;
for i = 1:length(d2.offM_highestsnr_fd01s)
    cellid = d2.offM_highestsnr_fd01s{i};
    if isempty(cellid), continue; end
    
    subplot(4,5,i);
    rasterphli(d2, cellid, d2.triggers(1:6:end));
    axis ij;
end
clear i cellid


%%
piece = '2012-09-18-3';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d1 = load_data(fullfile('/snle/acquisition/', piece, '/data005/data005'), loadopts);
d2 = load_data(fullfile('/snle/acquisition/', piece, '/data006/data006'), loadopts);

d1 = load_sta(d1, 'load_sta', []);
d1 = get_sta_summaries(d1, {4}, 'robust_std_method', 6);
d1 = calc_rf_snrs(d1);
[sorted, inds] = sortlownans(d1.stas.maxsigs);
d1.stas.offM_snrsorted = d1.cell_ids(inds(~isnan(sorted)));
clear sorted inds

cell_list_map = map_ei(d1, d2, 'master_cell_type', {4});
d2.offM_highestsnr_fd01s = cell_list_map(get_cell_indices(d1, d1.stas.offM_snrsorted(end:-1:end-20)));
clear cell_list_map;

% Check trigger spacing
plot(diff(d2.triggers), 'o-');

f = figure;
for i = 1:length(d2.offM_highestsnr_fd01s)
    cellid = d2.offM_highestsnr_fd01s{i};
    if isempty(cellid), continue; end
    
    subplot(4,5,i);
    rasterphli(d2, cellid, d2.triggers(1:6:end));
    axis ij;
    drawnow
end
clear i cellid