% Basics
piece = '2012-09-24-1';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'keep_vars'; 'd'};


%% Streamed runs and cell picks
streamedopts = loadopts;
streamedopts.load_sta = true;
streamedopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
streamedopts.set_polarities = {'guess', true};

d.d03s = load_data([piece '/data003'], streamedopts);
d.d03s.offM = [346  421 1232 1861 2807 2837 3766 4186 5207 5266 5761 6095 6826 6871 7397];

d.d05s = load_data([piece '/data005'], streamedopts);
d.d05s.offM = [346 1158 1951 2356 3121 4006 4126 5118 5266 5446 5761 6091 6676 7081 7397];
d.d05s.jfcones = {[]
    [112, 140, 106, 119]
    []
    [551, 531, 538, 556]
    [764, 779, 860, 859]
    []
    [1549, 1572, 1509, 1531]
    []
    [1700, 1709, 1717, 1691]
    []
    []
    []
    []
    []
    []};


%% C/R data006cm offM
crprintpath = [];%printpath('cr', piece);

if ~isfield(d, 'd06cm'), d.d06cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data006/data006'], loadopts); end
if ~isfield(d, 'd07cm'), d.d07cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data007/data007'], streamedopts); end

if ~isfield(d.d06cm, 'offM_fd03s')
    if ~isfield(d, 'd03_05_06_07_08_09_10'), d.d03_05_06_07_08_09_10 = load_data([piece '/d03-05-06-07-08-09-10-norefit'], loadopts); end
    cell_list_map = map_ei(d.d03s, d.d03_05_06_07_08_09_10, 'master_cell_type', {4});
    d.d06cm.offM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.offM));
    clear cell_list_map;
end

conerun = d.d03s;
conerun.rgcs = conerun.offM;
rasterrun = d.d06cm;
rasterrun.rgcs = rasterrun.offM_fd03s;
stablerun = d.d07cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d06cm d.d03s d.d07cm] = deal(cleared{:});
% 3x3 conerun just doesn't look that nice...


%% U/D data006cm offM
udprintpath = [];%printpath('ud', piece);

d.d03_05_06_07_08_09_10 = load_data([piece '/d03-05-06-07-08-09-10-norefit/d03-05-06-07-08-09-10-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data007/data007'], loadopts);

cell_list_map = map_ei(d.d03s, d.d03_05_06_07_08_09_10, 'master_cell_type', 'all');
d.d03_05_06_07_08_09_10.offM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.offM));
d.d06cm.offM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.offM));
d.d07cm.offM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.offM));
clear cell_list_map;
 
d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);
 
conerun = d.d03s;
conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
conerun = set_polarities(conerun);
conerun = get_sta_fits_from_vision(conerun);
conerun = load_cones(conerun, 'acquisition');
conerun = make_mosaic_struct(conerun);
conerun = make_voronoi_masks(conerun);
conerun.wwrgcs = conerun.offM;
 
stablerun = d.d07cm;
stablerun = load_sta(stablerun, struct('load_sta', [], 'guess_stimulus', false));
stablerun = set_polarities(stablerun, 'cell_specs', {{1},{8}});
stablerun = get_sta_fits_from_vision(stablerun);
stablerun.wwrgcs = stablerun.offM_fd03s;
 
ud06cm.wwrun = d.d06cm;
ud06cm.wwrun.wwrgcs = ud06cm.wwrun.offM_fd03s;
ud06cm.conerun = conerun;
ud06cm.stablerun = stablerun;
ww_compound_plot(ud06cm, 'all', 'triggers', ud06cm.wwrun.triggers(1:2:end), 'printpath', udprintpath);
% clear conerun stablerun udprintpath


% For J Freeman subunits analysis
save_udraster_data(ud06cm, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


%% C/R data008cm offM
crprintpath = printpath('cr', piece);

if ~isfield(d, 'd08cm'), d.d08cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data008/data008'], loadopts); end
if ~isfield(d, 'd10cm'), d.d10cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data010/data010'], streamedopts); end

if ~isfield(d.d08cm, 'offM_fd05s')
    if ~isfield(d, 'd03_05_06_07_08_09_10'), d.d03_05_06_07_08_09_10 = load_data([piece '/d03-05-06-07-08-09-10-norefit'], loadopts); end
    cell_list_map = map_ei(d.d05s, d.d03_05_06_07_08_09_10, 'master_cell_type', {4});
    d.d08cm.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
    clear cell_list_map;
end

conerun = d.d05s;
conerun.rgcs = conerun.offM;
rasterrun = d.d08cm;
rasterrun.rgcs = rasterrun.offM_fd05s;
stablerun = d.d10cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d08cm d.d05s d.d10cm] = deal(cleared{:});
% Not super clean, but maybe usable for BIG OFF midgets


%% U/D data008cm offM
udprintpath = [];%printpath('ud', piece);

d.d03_05_06_07_08_09_10 = load_data([piece '/d03-05-06-07-08-09-10-norefit/d03-05-06-07-08-09-10-norefit'], loadopts);
d.d08cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data008/data008'], loadopts);
d.d10cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data010/data010'], loadopts);

cell_list_map = map_ei(d.d05s, d.d03_05_06_07_08_09_10, 'master_cell_type', 'all');
d.d03_05_06_07_08_09_10.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d08cm.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d10cm.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
clear cell_list_map;
 
d.d08cm = read_stim_lisp_output(d.d08cm);
d.d08cm.stimulus = parse_stim_rgbs(d.d08cm.stimulus);
 
conerun = d.d05s;
conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
conerun = set_polarities(conerun);
conerun = get_sta_fits_from_vision(conerun);
conerun = load_cones(conerun, 'acquisition');
conerun = make_mosaic_struct(conerun);
conerun = make_voronoi_masks(conerun);
conerun.wwrgcs = conerun.offM;
 
stablerun = d.d10cm;
stablerun = load_sta(stablerun, struct('load_sta', [], 'guess_stimulus', false));
stablerun = set_polarities(stablerun, 'cell_specs', {{1}, {8}});
stablerun = get_sta_fits_from_vision(stablerun);
stablerun.wwrgcs = stablerun.offM_fd05s;
 
ud08cm.wwrun = d.d08cm;
ud08cm.wwrun.wwrgcs = ud08cm.wwrun.offM_fd05s;
ud08cm.conerun = conerun;
ud08cm.stablerun = stablerun;
ww_compound_plot(ud08cm, 'all', 'triggers', ud08cm.wwrun.triggers(1:2:end), 'printpath', udprintpath);
% clear conerun stablerun udprintpath


% For J Freeman subunits analysis
save_udraster_data(ud08cm, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


%% U/D data009cm offM
udprintpath = printpath('ud', piece);

d.d03_05_06_07_08_09_10 = load_data([piece '/d03-05-06-07-08-09-10-norefit/d03-05-06-07-08-09-10-norefit'], loadopts);
d.d09cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data009/data009'], loadopts);
d.d10cm = load_data([piece '/d03-05-06-07-08-09-10-norefit/data010/data010'], loadopts);

cell_list_map = map_ei(d.d05s, d.d03_05_06_07_08_09_10, 'master_cell_type', 'all');
d.d03_05_06_07_08_09_10.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d09cm.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d10cm.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
clear cell_list_map;
 
d.d09cm = read_stim_lisp_output(d.d09cm);
d.d09cm.stimulus = parse_stim_rgbs(d.d09cm.stimulus);
 
conerun = d.d05s;
conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
conerun = set_polarities(conerun);
conerun = get_sta_fits_from_vision(conerun);
conerun = load_cones(conerun, 'acquisition');
conerun = make_mosaic_struct(conerun);
conerun = make_voronoi_masks(conerun);
conerun.wwrgcs = conerun.offM;
 
stablerun = d.d10cm;
stablerun = load_sta(stablerun, struct('load_sta', [], 'guess_stimulus', false));
stablerun = set_polarities(stablerun, 'cell_specs', {{1}, {8}});
stablerun = get_sta_fits_from_vision(stablerun);
stablerun.wwrgcs = stablerun.offM_fd05s;
 
ud09cm.wwrun = d.d09cm;
ud09cm.wwrun.wwrgcs = ud09cm.wwrun.offM_fd05s;
ud09cm.conerun = conerun;
ud09cm.stablerun = stablerun;
ww_compound_plot(ud09cm, 'all', 'triggers', ud09cm.wwrun.triggers(1:2:end), 'printpath', udprintpath);
% clear conerun stablerun udprintpath


% For J Freeman subunits analysis
save_udraster_data(ud09cm, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));
