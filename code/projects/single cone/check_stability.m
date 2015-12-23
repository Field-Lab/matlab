%% Stability check / Show cone maps for RGC selection
piece = '/snle/acquisition/2013-03-11-1/';
runA = 'data002';
runB = 'data018';
conesA = {};
conesB = {};

datarunA = load_data(fullfile(piece, runA, runA));
datarunA = load_sta(datarunA, 'load_sta', []);
datarunA = load_params(datarunA);
datarunA = load_cones(datarunA, conesA{:});
datarunA = make_mosaic_struct(datarunA);
datarunA = get_sta_fits_from_vision(datarunA);

datarunB = load_data(fullfile(piece, runB, runB));
datarunB = load_sta(datarunB, 'load_sta', []);
datarunB = load_params(datarunB);
datarunB = load_cones(datarunB, conesB{:});
datarunB = make_mosaic_struct(datarunB);
datarunB = get_sta_fits_from_vision(datarunB);

% datarunC = load_data('/Volumes/War/Acquisition/2011-07-05-2/data004/data004');
% datarunC = load_sta(datarunC, 'load_sta', []);
% datarunC = load_params(datarunC);
% conesC = '_Volumes_War_Acquisition_2011-07-05-2_data004_data004-bayes-msf_15.00-BW-2-4';
% datarunC = import_single_cone_data(datarunC, conesC);
% datarunC = get_sta_fits_from_vision(datarunC);


% Export for overlays in Intaglio
% plot_cone_mosaic(datarunB, 'fig_or_axes', 1,...
%                 'cone_size', 6, 'cone_colors', [0 0 0])
% print(1, '~/Desktop/cone-mosaic-a.pdf','-dpdf')
%  
% plot_cone_mosaic(datarunA, 'fig_or_axes', 2,...
%                 'cone_size', 3, 'cone_colors', [1 1 1])
% print(2, '~/Desktop/cone-mosaic-b.pdf','-dpdf')


% Overlay in Matlab
overlay_cone_mosaics(datarunA, datarunB);
%overlay_cone_mosaics(datarunA, datarunB, 'scale1', [2 2], 'scale2', [3 3]);

% Show if any RGCs we picked were stable FUCKERS!!!!
%plot_rf_summaries(datarunA, nullf13, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')
