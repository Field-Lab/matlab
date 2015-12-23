%% Basics
piece = '2012-09-18-3';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'udprintpath'; 'allconesprintpath'; 'keep_vars'; 'd'};

%% Stability check for Jeremy Freeman's repeat analysis
d.d03s = load_data([piece '/streamed/data003/data003'], staopts);
d.d03s = load_cones(d.d03s);
d.d05s = load_data([piece '/streamed/data005/data005'], staopts);
d.d05s = load_cones(d.d05s);
overlay_cone_mosaics(d.d03s, d.d05s);
plot_rf_summaries(d.d03s, {4}, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y');