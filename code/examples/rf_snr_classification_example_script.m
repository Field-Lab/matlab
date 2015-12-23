%% Spitting out deciles to Vision classification.txt
datarun = load_data('/snle/acquisition/2011-05-11-6/data009/data009');
tic; datarun = load_sta(datarun, 'save_rf', true); toc
rf_snr_classification(datarun);


%% Getting the highest SNR cells from a class
piece = '2012-08-09-1';
data = 'streamed/data001/data001';
cellspec = {4};

robust_std_method = 6;
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

datarun = load_data(fullfile(piece, data), loadopts);
datarun = load_sta(datarun, 'load_sta', []);
datarun = get_sta_summaries(datarun, cellspec, 'robust_std_method', robust_std_method);
datarun = calc_rf_snrs(datarun);
[sorted, inds] = sortlownans(datarun.stas.medsigs);
snr_sorted_cells = datarun.cell_ids(inds(~isnan(sorted)));
clear sorted inds

highest_snr_cells = snr_sorted_cells(end:-1:end-10)'