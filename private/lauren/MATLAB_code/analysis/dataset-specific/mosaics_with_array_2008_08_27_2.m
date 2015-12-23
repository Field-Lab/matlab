clear all


datarun = load_data('2008-08-27-2/data001-lh/data001-lh');


datarun = load_index(datarun);
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));
datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
datarun = load_ei(datarun,[]);
datarun = get_sta_fits_from_vision(datarun,'all');


% compute the transformation using the clicked points in datarun.piece.array
datarun = compute_monitor_to_array_transformation(datarun);




plot_rf_summaries(datarun, {1}, 'fit_width', 1, 'array', true)
hold on
plot_rf_summaries(datarun, {3}, 'fit_width', 1)

%plot stim electrodes 46 and 49
plot([datarun.ei.position_sta(46, 1) datarun.ei.position_sta(49, 1)], [datarun.ei.position_sta(46, 2) datarun.ei.position_sta(49, 2)], 'k.')

%% other stuff from example file monitor_alignment.m

% plot an RF with the array
plot_rf(datarun,754,'array',true, 'foa', 0)
% add the EI
hold on;plot_ei(datarun,754,'coordinates','sta','pretty_axes',0,'alpha',0)

% plot an RF with the array
plot_rf(datarun,677,'array',true, 'foa', 0)
% add the EI
hold on;plot_ei(datarun,677,'coordinates','sta','pretty_axes',0,'alpha',0)

% plot an RF with the array
plot_rf(datarun,676,'array',true, 'foa', 0)
% add the EI
hold on;plot_ei(datarun,676,'coordinates','sta','pretty_axes',0,'alpha',0)