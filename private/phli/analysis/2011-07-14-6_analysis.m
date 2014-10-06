piece = '2011-07-14-6';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'udprintpath'; 'allconesprintpath'; 'keep_vars'; 'd'};


%% streamed/data001
off_midgets_d01s = [53 76 98 118 122 139 190 215 238];

d01s = load_data('2011-07-14-6/streamed/data001/data001');
d01s = load_params(d01s);
d01s = load_ei(d01s, {4});

d01s = load_sta(d01s,'load_sta',[]);
d01s = set_polarities(d01s);
d01s = get_sta_fits_from_vision(d01s);

d01s = load_cones(d01s);
d01s = make_mosaic_struct(d01s);
d01s = make_voronoi_masks(d01s);


%% streamed/data005
d.d05s = load_data('2011-07-14-6/streamed/data005/data005', loadopts);
d.d05s = load_sta(d.d05s, 'load_sta', []);
d.d05s = set_polarities(d.d05s);
d.d05s = load_cones(d.d05s, 'acquisition');
d.d05s = make_mosaic_struct(d.d05s);
d.d05s = make_voronoi_masks(d.d05s);
d.d05s = get_sta_fits_from_vision(d.d05s);

% Used for allcone run, but probably data are no good.
d.d05s.offM = [68 169 199 203 296];


%% data004
d04 = load_data('2011-07-14-6/data004');
d04 = load_neurons(d04);
d04 = load_ei(d04, 'all');

d04.triggers = d04.triggers(1:3:end);


%% data005
d05 = load_data('2011-07-14-6/data005/data005');
d05 = load_params(d05);
d05 = load_sta(d05, 'load_sta', []);
d05 = load_cones(d05, '^2');
d05 = make_mosaic_struct(d05);
d05 = load_ei(d05, 'all');


%% data006
d06 = load_data('2011-07-14-6/data006');
d06 = load_neurons(d06);
d06 = load_ei(d06, 'all');

d06.triggers = d06.triggers(1:3:end);


%% data009
d09 = load_data('2011-07-14-6/data009');
d09 = load_neurons(d09);
d09 = load_ei(d09, 'all');


%% stimuli
s04 = read_stim_lisp_output('2011-07-14-6/s04');
s06 = read_stim_lisp_output('2011-07-14-6/s06');
d09 = read_stim_lisp_output(d09, ':2011-07-14-6:2011-07-14-6_f05_allcones');

s04.mapdir = [server_data_path '2011-07-14-6/2011-07-14-6_f01_1and2'];
s06.mapdir = [server_data_path '2011-07-14-6/2011-07-14-6_f01_12serial'];

mapstimuli = {s04 s06};
for s = 1:length(mapstimuli)
    stim = mapstimuli{s};
    mapims = {};
    mapnyc = {};
    for i = stim.mapindices
        fprintf('.');
        mapims{i+1} = sparse(dlmread([stim.mapdir sprintf('/map-%.4d.txt', i)]));
        mapnyc{i+1} = map2manhattan(mapims{i+1}, 'quiet', true);
    end; 
    mapstimuli{s}.mapims = mapims;
    mapstimuli{s}.mapnyc = mapnyc;
    fprintf('\n');
end; clear s i stim mapims mapnyc;

s04 = mapstimuli{1};
s06 = mapstimuli{2};
s09 = mapstimuli{3};
clear mapstimuli;


%% Map d04 to d01s on EIs, get handpicked off midgets
cell_list_map = map_ei(d01s, d04, 'master_cell_type', {4});
off_midgets_fd01s_d04 = cell_list_map(get_cell_indices(d01s, off_midgets_d01s));
clear cell_list_map;


%% Map d06 to d01s on EIs, get handpicked off midgets
cell_list_map = map_ei(d01s, d06, 'master_cell_type', {4});
off_midgets_fd01s_d06 = cell_list_map(get_cell_indices(d01s, off_midgets_d01s));
clear cell_list_map;


%% Map d05 to d05s on EIs, get handpicked off midgets
cell_list_map = map_ei(d05s, d05, 'master_cell_type', {4});
off_midgets_fd05s_d05 = cell_list_map(get_cell_indices(d05s, off_midgets_d05s));
clear cell_list_map;

% Handcleaned... No 68...
off_midgets_fd05s_d05 = {[] 4951 6511 6544 3033};


%% Map d09 to d05s on EIs, get handpicked off midgets
cell_list_map = map_ei(d05s, d09, 'master_cell_type', {4});
d09.offM_fd05s = cell_list_map(get_cell_indices(d05s, d05s.offM));
clear cell_list_map;

% The upshot, matching d05s cell id 169
off_midgets_fd05s_d09 = {[] 4773 [] [] []};


%% d09-fd05 doesn't seem to work.
d09fd05 = load_data('2011-07-14-6/data009-from-data005', loadopts);
cell_list_map = map_ei(d05s, d09fd05, 'master_cell_type', {4}, 'electrode_threshold', 0, 'corr_threshold', 0.5);
d09fd05.offM_fd05s = cell_list_map(get_cell_indices(d05s, d05s.offM));


%% Build additivity structs
additivity_d04.conerun = d01s;
additivity_d04.conerun.addrgcs = off_midgets_d01s;
additivity_d04.addrun  = d04;
additivity_d04.addstim = s04;
additivity_d04.addrun.addrgcs  = off_midgets_fd01s_d04;
additivity_d04.stablerun = d05;


%% Plot rasters d04: superposition, contrast-response
additivity_compound_plot(additivity_d04, 'all');
% index 2 == 76 == 1112,  stability not clear
% index 8 == 215 == 2646, no stability data
% others look pretty good


%% Plot rasters d06: crosstalk control
crosstalk_raster_plot(d06, s06, off_midgets_fd01s_d06);


%% Contrast response analysis
cr04 = additivity_contrast_response({d04 d06}, {s04 s06}, {off_midgets_fd01s_d04 off_midgets_fd01s_d06}, 0.75, 0.02, 'start', -0.1);
cr04(:,4,:) = cr04(:,1,:) + cr04(:,2,:);


%% Compare C/R scaled by STA weights to least squares
crstruct04.crfs = cr04;
crstruct04.crrun = d04;
crstruct04.crstim = s04;
crstruct04.crcells = off_midgets_fd01s_d04;
crstruct04.conerun = d01s;
crstruct04.conecells = off_midgets_d01s;
crstruct04.stablerun = d05;
inds = {[2 4] [6 8] [9]};
for i = 1:length(inds)
    figure;
    cr_compound_plot(crstruct04, inds{i})
end; clear i inds


%% Allcones data009
% Sucks; maybe double check that we're getting the right RGCs and cones?
allconesprintpath = [];%printpath('allcones', piece);

if ~exist('d') || ~isfield(d, 'd05cm')
    d.d05cm = load_data([piece '/d05-07-08-09-10-n/d05-f-05-07-08-09-10/d05-f-05-07-08-09-10'], loadopts);
end
if ~isfield(d, 'd09cm')
    d.d09cm = load_data([piece '/d05-07-08-09-10-n/d09-f-05-07-08-09-10/d09-f-05-07-08-09-10'], loadopts);
end
cell_list_map = map_ei(d.d05s, d.d05cm, 'master_cell_type', {4});
d.d09cm.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
clear cell_list_map;
d.d09cm.offM_fd05s{5} = 3033; % Found by setting corr threshold down to 0.75; the cell in d05s is contaminated but looks clean on d05cm!

d.d09cm = read_stim_lisp_output(d.d09cm, '2011-07-14-6_f05_allcones');
d.d09cm.stimulus = parse_stim_rgbs(d.d09cm.stimulus);

conerun = d.d05s;
datarun = d.d09cm;
conerun.rgcs = conerun.offM;
datarun.rgcs = datarun.offM_fd05s;
datarun.stimulus.triggers = datarun.triggers(1:3:end);

conerun = load_sta(conerun, 'load_sta', []);
conerun = get_sta_fits_from_vision(conerun);
conerun = set_polarities(conerun);

if ~isfield(d, 'd10cm')
    d.d10cm = load_data([piece '/d05-07-08-09-10-n/d10-f-05-07-08-09-10/d10-f-05-07-08-09-10'], loadopts);
end
stablerun = d.d10cm;
stablerun.rgcs = d.d09cm.offM_fd05s;
stablerun = load_sta(stablerun, 'load_sta', []);
stablerun = set_polarities(stablerun, 'cell_specs', {{1} {8}});
allcones_plot(datarun, conerun, 'printpath', allconesprintpath, 'urgbs', datarun.stimulus.urgbs(1), 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});


%% Scatter plot DEPRECATED
% figure; hold on;
% plot(squeeze(cr04(:,4,:)), squeeze(cr04(:,3,:)), 'b.');
% plot(-10:320, -10:320, 'k--');
% axis equal tight;
% title 'Off Midgets 2011-07-14-6';
% xlabel 'Linear Prediction';
% ylabel 'Observed Firing Rate';


%% Scatter plot added to 2011-06-24-6 DEPRECATED
% subplot(1, 2, 1);
% h201107146 = plot(squeeze(cr04(:,4,:)), squeeze(cr04(:,3,:)), 'r.');
% legend([h1(1) h2(1)], '2011-06-24-6', '2011-07-14-6');


%% Attempt to look at one cell from all-cones run
% d09_fixtrigg = d09;
% d09_fixtrigg.triggers = d09.triggers(1:3:end);
% crosstalk_raster_plot(d09_fixtrigg, d09.stimulus, d09.offM_fd05s, 'maps', rot90([0:9;10:19;20:29],2));