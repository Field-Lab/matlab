% Basics
piece = '2012-09-06-0';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%% Streamed runs and cell picks


d.d04s = load_data([piece, '/streamed/data004/data004'], staopts);
d.d04s.offMjf = [2146 2491 2812 3766 4581 4936 5103 5746 6004 6200 6650 7008];

% Processed as in the new SINGLE_CONE_VORONOI_MAPS_SCRIPT_PHLI
d.d04s.offMjfcones = [
    302, 322, 356, 366
    654, 674, 673, 693
    531, 565, 556, 601
    967, 1007, 1032, 1014
    994, 1015, 1016, 995
    1144, 1159, 1143, 1192
    873, 895, 910, 896
    1115, 1116, 1122, 1133
    1000,1001,1019, 1036
    882, 898, 906, 942
    797, 806, 832, 833
    636, 645, 649, 627
];

% Used for allcones?
d.d04s.offMallcones = [275 436 904 932 1892 2146 2716 2812 2986 3212 3395 3722 4126 4936 5103 5746 6031 6289 6317 6889 6947 7022 7505];


%% C/R/U/D data007cm
crprintpath = [];%printpath('cr', piece);

if ~isfield(d, 'd04cm'), d.d04cm = load_data([piece '/d04-07-09-norefit/data004/data004'], staopts); end
if ~isfield(d, 'd07cm'), d.d07cm = load_data([piece '/d04-07-09-norefit/data007/data007'], loadopts); end
if ~isfield(d, 'd08cm'), d.d09cm = load_data([piece '/d04-07-09-norefit/data009/data009'], staopts); end

if ~isfield(d.d07cm, 'offM_fd04s')
    if ~isfield(d, 'd04_07_09'), d.d04_07_09 = load_data([piece '/d04-07-09-norefit'], loadopts); end
    cell_list_map = map_ei(d.d04s, d.d04_07_09);
    d.d07cm.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offMjf));
    clear cell_list_map;
end

rasterrun = d.d07cm;
rasterrun.rgcs = rasterrun.offM_fd04s;
conerun = d.d04cm;
conerun.rgcs = rasterrun.rgcs;
stablerun = d.d09cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d07cm d.d04cm d.d09cm] = deal(cleared{:});



udprintpath = [];%printpath('ud', piece);
rasterrun.wwrgcs = rasterrun.offM_fd04s;
conerun.wwrgcs = rasterrun.wwrgcs;
stablerun.wwrgcs = rasterrun.wwrgcs;
ud07cm.wwrun = rasterrun;
ud07cm.conerun = conerun;
ud07cm.stablerun = stablerun;
ww_compound_plot(ud07cm, 'all', 'triggers', ud07cm.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear conerun stablerun udprintpath



%% For J Freeman subunits analysis
if ~isfield(d, 'd07'), d.d07 = load_data([piece '/data007'], loadopts); end
if ~isfield(d.d07, 'offM_fd04s')
    cell_list_map = map_ei(d.d04s, d.d07);
    d.d07.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offMjf));
    clear cell_list_map;
end
if ~isfield(d, 'd09'), d.d09 = load_data([piece '/data009'], staopts); end
if ~isfield(d.d09, 'offM_fd04s')
    cell_list_map = map_ei(d.d04s, d.d09);
    d.d09.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offMjf));
    clear cell_list_map;
end

d.d07 = read_stim_lisp_output(d.d07);
d.d07.stimulus = parse_stim_rgbs(d.d07.stimulus);
ud07.wwrun = d.d07;
ud07.wwrun.wwrgcs = ud07.wwrun.offM_fd04s;
ud07.conerun = d.d04s;
ud07.conerun.wwrgcs = ud07.conerun.offMjf;
ud07.stablerun = d.d09;
ud07.stablerun.wwrgcs = d.d09.offM_fd04s;
save_udraster_data(ud07, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


%% data008 allconesd04, can also check concatenated
allconesprintpath = [];%printpath('allcones', piece);

if ~isfield(d, 'd08'), d.d08 = load_data([piece '/data008'], loadopts); end
d.d08 = read_stim_lisp_output(d.d08);

cell_list_map = map_ei(d.d04s, d.d08, 'master_cell_type', 'all');
d.d08.offM_f04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offMallcones));
clear cell_list_map;

conerun = d.d04s;
datarun = d.d08;
conerun.rgcs = conerun.offMallcones;
datarun.rgcs = datarun.offM_f04s;
datarun.stimulus.triggers = datarun.triggers(1:2:end);

if ~isfield(d, 'd09'), d.d09 = load_data([piece '/data009'], staopts); end
cell_list_map = map_ei(d.d04s, d.d09, 'master_cell_type', 'all');
d.d09.offM_f04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offMallcones));
clear cell_list_map;

stablerun = d.d09;
stablerun.rgcs = stablerun.offM_f04s;

allcones_plot(datarun, conerun, 'coneind', 'localmax', 'urgbs', {[1 1 1].*-0.288 [1 1 1].*-0.48}, 'printpath', allconesprintpath, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});

leave(keep_vars);