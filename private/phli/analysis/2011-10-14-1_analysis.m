%% Basic load data
piece = '2011-10-14-1';
udprintpath = [printpath '/ud'];

d01s = load_data('2011-10-14-1/streamed/data001-0/data001-0');
d03 = load_data('2011-10-14-1/data003');
d04 = load_data('2011-10-14-1/data004');
d06 = load_data('2011-10-14-1/data006');

dataruns = {'d01s' 'd03' 'd04' 'd06'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% EIs
dataruns = {'d01s' 'd03' 'd04' 'd06'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_ei(datarun, 'all');
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i



%% STAs and cones
dataruns = {'d01s' 'd04' 'd06'};
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
offmidgets_d01s = [136 947 2311 7471 3063 3797 4996 5717 6811];


%% Stimuli
d03 = read_stim_lisp_output(d03);

% % Wiggidy-wack
% s03 = read_stim_lisp_output(piece, 's03');
% s05 = read_stim_lisp_output(piece, 's05');
% 
% 
% stimuli = {'s03' 's05'};
% for i = 1:length(stimuli)
%     assignin('caller', stimuli{i}, stim);
% end; clear stim stimuli i mapdirs


%% Map d03 to d01s on EIs, get handpicked off midgets
cell_list_map = map_ei(d01s, d03, 'master_cell_type', {4});
offmidgets_fd01s_d03 = cell_list_map(get_cell_indices(d01s, offmidgets_d01s));
clear cell_list_map;


%% Map d04 to d01s on EIs, get handpicked off midgets
cell_list_map = map_ei(d01s, d04, 'master_cell_type', {4});
offmidgets_fd01s_d04 = cell_list_map(get_cell_indices(d01s, offmidgets_d01s));
clear cell_list_map;

% the upshot
offmidgets_fd01s_d04 = {[7741] [946] [2298] [7547] [3068] [3781] [5056] [5716] [6811]};


%% Map d06 to d01s on EIs, get handpicked off midgets
cell_list_map = map_ei(d01s, d06, 'master_cell_type', {4});
offmidgets_fd01s_d06 = cell_list_map(get_cell_indices(d01s, offmidgets_d01s));
clear cell_list_map;

% the upshot
offmidgets_fd01s_d06 = {[31] [946] [2206] [7471] [] [] [5056] [5716] [6811]};


%%
ww03.wwrun = d03;
ww03.wwrun.wwrgcs = offmidgets_fd01s_d03;
ww03.conerun = d01s;
ww03.conerun.wwrgcs = offmidgets_d01s;
ww03.stablerun = d04;
ww03.stablerun.wwrgcs = offmidgets_fd01s_d04;
ww_compound_plot(ww03, 'all', 'triggers', d03.triggers(1:2:end), 'printpath', udprintpath);

% 1, 3 are a little interesting...

% For J Freeman subunits analysis
save_udraster_data(ww03, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


%% For figures
ww_compound_plot(ww03, 9, 'triggers', d03.triggers(1:2:end), 'scaled_up', 10, 'MarkerSize', 15, 'set_ylims', [0 55]);
ww_compound_plot(ww03, 4, 'triggers', d03.triggers(1:2:end), 'scaled_up', 10, 'MarkerSize', 15, 'set_ylims', [0 45], 'az_pad_factor', 3);