% Basics
piece = '2012-08-21-2';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'udprintpath'; 'allconesprintpath'; 'keep_vars'; 'd'};

%%
d.d01s = load_data([piece '/streamed/data001/data001'], loadopts);
d.d01s.offM = [153 196 391 408 887 1268 1777 2223 2266 2708 2828 3033 4956 5447 6005 6721 7682];

%% data003 Allcones analysis
allconesprintpath = printpath('allcones');

if ~exist('d') || ~isfield(d, 'd03')
    d.d03 = load_data([piece '/data003'], loadopts);
end
d.d03 = read_stim_lisp_output(d.d03);

cell_list_map = map_ei(d.d01s, d.d03, 'master_cell_type', 'all');
d.d03.offM_f01s = cell_list_map(get_cell_indices(d.d01s, [d.d01s.offM]));
clear cell_list_map;

conerun = d.d01s;
datarun = d.d03;
conerun.rgcs = d.d01s.offM;
datarun.rgcs = d.d03.offM_f01s;
datarun.stimulus.triggers = datarun.triggers(1:2:end);

conerun = load_sta(conerun, 'load_sta', []);
conerun = get_sta_fits_from_vision(conerun);
conerun = set_polarities(conerun);

if ~isfield(d, 'd04')
    d.d04 = load_data([piece '/data004'], loadopts);
end
cell_list_map = map_ei(d.d01s, d.d04, 'master_cell_type', 'all');
d.d04.offM_f01s = cell_list_map(get_cell_indices(d.d01s, [d.d01s.offM]));
clear cell_list_map;

stablerun = d.d04;
stablerun.rgcs = d.d04.offM_f01s;
stablerun = load_sta(stablerun, 'load_sta', []);
stablerun = set_polarities(stablerun, 'cell_specs', {{1} {8}});

allcones_plot(datarun, conerun, 'coneind', 'localmax', 'urgbs', {[1 1 1].*-0.36 [1 1 1].*-0.48}, 'printpath', allconesprintpath, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});


leave(keep_vars);