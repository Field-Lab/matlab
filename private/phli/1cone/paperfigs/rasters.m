% Started rewriting before found simplerasters.m; probably go back to
% simplerasters.m

% Basics
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%%
piece = '2012-09-13-2';
d.d05s = load_data(fullfile(piece, '/streamed/data005/data005'), staopts);
d.d05 = load_data([piece '/data005'], staopts);
d.d07 = load_data([piece '/data007'], loadopts);
d.d07 = read_stim_lisp_output(d.d07);
d.d07.stimulus = parse_stim_rgbs(d.d07.stimulus);
d.d07.cellmaps = containers.Map('KeyType', 'char', 'ValueType', 'any');
d.d07.cellmaps('d05') = map_ei(d.d05, d.d07);

figure;
sanesubplot(3, 2, {1:2 1:2});
plot_rf_stimmap(d.d05, 2116, d.d07.stimulus.mapnycpoly{1}([3 1]), 'fit', false, 'colors', ['kkkk']', 'scaled_up', 10);
axis([88.5 116.5 82.5 110.5])
urgbs = num2cell(find(d.d07.stimulus.urgb.singles & sum(d.d07.stimulus.urgb.intensities) == -1.44));
urgbs = {[] []; [] []; urgbs{[2 4]}};
urgb_raster_subplot(d.d07, 2117, d.d07.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k');


figure;
sanesubplot(3, 2, {1:2 1:2});
plot_rf_stimmap(d.d05, 6992, d.d07.stimulus.mapnycpoly{1}([4 2]), 'fit', false, 'colors', ['kkkk']', 'scaled_up', 10);
urgbs = num2cell(find(d.d07.stimulus.urgb.singles & sum(d.d07.stimulus.urgb.intensities) == -1.44));
urgbs = {[] []; [] []; urgbs{[1 3]}};
urgb_raster_subplot(d.d07, 6995, d.d07.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k');


figure;
sanesubplot(3, 2, {1:2 1:2});
plot_rf_stimmap(d.d05, 4441, d.d07.stimulus.mapnycpoly{1}([3 1]), 'fit', false, 'colors', ['kkkk']', 'scaled_up', 10);
x = 215; y = 101; del = 70; axis([x x+del y y+del]);
urgbs = num2cell(find(d.d07.stimulus.urgb.singles & sum(d.d07.stimulus.urgb.intensities) == -1.44));
urgbs = {[] []; [] []; urgbs{[2 4]}};
urgb_raster_subplot(d.d07, 4517, d.d07.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k');


%%
piece = '2011-12-13-2';
d.d04s = load_data(fullfile(piece, 'streamed', 'data004-0', 'data004-0'), staopts);
d.d08s = load_data(fullfile(piece, 'streamed', 'data008-0', 'data008-0'), staopts);
d.d11s = load_data(fullfile(piece, 'streamed', 'data011-0', 'data011-0'), staopts);

d.d12 = load_data(fullfile(piece, 'data012'), loadopts);
d.d12 = read_stim_lisp_output(d.d12);
d.d12.stimulus = parse_stim_rgbs(d.d12.stimulus);
d.d12.map08s = map_ei(d.d08s, d.d12);
d.d12.map11s = map_ei(d.d11s, d.d12);

figure;
sanesubplot(3, 2, {1:2 1:2});
conestoplot = [4 3];
d.d08s.onMrf = get_rf(d.d08s, 856);
d.d11s.onMrf = get_rf(d.d11s, 857);
d.compound = d.d08s;
d.compound.names.short_name = '2011-12-13-2_compound_d08s+d11s';
d.compound.stas.rfs{get_cell_indices(d.compound, 856)} = d.d11s.onMrf;
plot_rf_stimmap(d.compound, 856, d.d12.stimulus.mapnycpoly{1}(conestoplot), 'fit', false, 'colors', ['kk']', 'scaled_up', 10);
x = 118.5; y = 183.5; del = 30; axis([x x+del y y+del]);
urgbs = num2cell(find(d.d12.stimulus.urgb.singles & sum(d.d12.stimulus.urgb.intensities(conestoplot,:)) == 1.44));
temp = cell(3, length(urgbs));
temp(3,:) = urgbs;
urgbs = temp;
urgb_raster_subplot(d.d12, 727, d.d12.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k');