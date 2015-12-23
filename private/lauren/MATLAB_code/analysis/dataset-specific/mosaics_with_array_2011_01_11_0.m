clear all


datarun = load_data('2011-01-11-0/data030-lh/data030-lh');


datarun = load_index(datarun);
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));
datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
datarun = load_ei(datarun,[]);
datarun = get_sta_fits_from_vision(datarun,'all');


% compute the transformation using the clicked points in datarun.piece.array
datarun = compute_monitor_to_array_transformation(datarun);

lineColors(1,:) = [90 156 0]/255; %pale grass
lineColors(2,:) = [255 124 59]/255; %salmon
lineColors(3,:) = [101 52 255]/255; %purple
lineColors(4,:) = [52 198 247]/255; %aqua
lineColors(5,:) = [238 55 128]/255; %calm magenta

%% plot full mosaics


plot_rf_summaries(datarun, {1}, 'fit_width', 1, 'array', true, 'fit_color', lineColors(1,:))
hold on
plot_rf_summaries(datarun, {2}, 'fit_width', 1, 'fit_color', lineColors(2,:))
plot_rf_summaries(datarun, {3}, 'fit_width', 1, 'fit_color', lineColors(3,:))
plot_rf_summaries(datarun, {4}, 'fit_width', 1, 'fit_color', lineColors(4,:))


%% plot individual RFs


plot_rf_summaries(datarun, [799], 'fit_width', 2, 'fit_color', lineColors(4,:))






%% other stuff from example file monitor_alignment.m
% 

%plot stim electrodes 46 and 49
%plot([datarun.ei.position_sta(46, 1) datarun.ei.position_sta(49, 1)],
%[datarun.ei.position_sta(46, 2) datarun.ei.position_sta(49, 2)], 'k.')


% % plot an RF with the array
% plot_rf(datarun,754,'array',true, 'foa', 0)
% % add the EI
% hold on;plot_ei(datarun,754,'coordinates','sta','pretty_axes',0,'alpha',0)
% 
% % plot an RF with the array
% plot_rf(datarun,677,'array',true, 'foa', 0)
% % add the EI
% hold on;plot_ei(datarun,677,'coordinates','sta','pretty_axes',0,'alpha',0)
% 
% % plot an RF with the array
% plot_rf(datarun,676,'array',true, 'foa', 0)
% % add the EI
% hold on;plot_ei(datarun,676,'coordinates','sta','pretty_axes',0,'alpha',0)