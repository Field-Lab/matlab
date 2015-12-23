%% Basic load data
piece = '2011-10-25-9';
udprintpath = [printpath() '/ud'];

d06s = load_data([piece '/streamed/data006-0/data006-0']);
d12s = load_data([piece '/streamed/data012-0/data012-0']);
d11 = load_data([piece '/data011']);

dataruns = {'d06s' 'd12s' 'd11'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% EIs
dataruns = {'d06s' 'd12s' 'd11'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_ei(datarun, 'all');
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% STAs and cones
dataruns = {'d06s' 'd12s'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = set_polarities(datarun);
    datarun = get_sta_fits_from_vision(datarun);
    
    datarun = load_cones(datarun, 1);
    datarun = make_mosaic_struct(datarun);
    datarun.cones.mosaic.voronoi_masks_scale = 3;
    datarun = make_voronoi_masks(datarun);
    
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% Picks
midgets_d06s = [181 619 1411 1966 1982 3422 4066 4489 1187 7760 6814 6093 4866 1096 5436 3258 2193];
% 3258 2193 are ON, rest OFF

%% Stimuli
d11 = read_stim_lisp_output(d11);


%% Map d11 to d06s on EIs, get handpicked on/off midgets
cell_list_map = map_ei(d06s, d11, 'master_cell_type', {3, 4});
midgets_fd06s_d11 = cell_list_map(get_cell_indices(d06s, midgets_d06s));
clear cell_list_map;


%% Map d12s to d06s on EIs, get handpicked on/off midgets
cell_list_map = map_ei(d06s, d12s, 'master_cell_type', {3, 4});
midgets_fd06s_d12s = cell_list_map(get_cell_indices(d06s, midgets_d06s));
clear cell_list_map;


%%

ww11.wwrun = d11;
ww11.wwrun.wwrgcs = midgets_fd06s_d11;
ww11.conerun = d06s;
ww11.conerun.wwrgcs = midgets_d06s;
ww11.stablerun = d12s;
ww11.stablerun.wwrgcs = midgets_fd06s_d12s;
ww_compound_plot(ww11, 'all', 'triggers', d11.triggers(1:2:end), 'printpath', udprintpath);

% Not much insight into ONs

% fig 13 (cell 15, id 5436) pretty nice
% fig 8 (cell 8, id 4489) nice
% fig 6 (cell 6, id 3422) pretty good, but stability is an issue

% fig 1 (cell 1, id 181) a little concerning the asymmetry and movement suggestion


% For J Freeman subunits analysis
save_udraster_data(ww11, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));



%% Separate plotting for figures
figure;
plot_voronoi_over_rf(ww11.conerun, ww11.conerun.wwrgcs(6), 'az_aspect_ratio', 1, 'scaled_up', 10);
scalemap = 1 / ww11.conerun.cones.mosaic.voronoi_masks_scale;
patchmanhattan(ww11.wwstim.mapnyc{1}, 'scale', [scalemap scalemap], 'colors', [1 0 0; 0 0.5 0; 0 0 1; 1 0 1], 'patchopts', {'LineWidth', 2});
clear scalemap;