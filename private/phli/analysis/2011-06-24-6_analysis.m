piece = '2011-06-24-6';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
udprintpath = [];%[printpath() '/ud'];
allconesprintpath = [];%[printpath() '/allcones'];
keep_vars = {'piece'; 'loadopts'; 'udprintpath'; 'allconesprintpath'; 'keep_vars'; 'd'};

% Stability 3-5, 8-13 decent
% Stability 5-8 good
% Stability 13-15 totally unusably terrible

%% Runs with picked cells
d.d05s = load_data([piece '/streamed/data005/data005'], staopts);
d.d05s = load_cones(d.d05s, 'cquisition');
d.d05s = make_voronoi_masks(d.d05s);
d.d05s.offM = [50 78 64 137 158 120 15 29 42];

d.d08s = load_data([piece '/streamed/data008/data008'], staopts);
d.d08s = load_cones(d.d08s, 'cquisition');
d.d08s = make_voronoi_masks(d.d08s);
d.d08s.offP = [34 75 88 96 125 198];
d.d08s.offM = [16 40 57 83 144 183]; % Allcones



%% Allcones d14, allcones_plot_old
% Stability is terrible; totally unusable!
%
% Useful for demonstrating movement issue?
% Didn't really pan out.

 
% d.d05 = load_data([piece '/data005'], loadopts);
% cell_list_map = map_ei(d.d05s, d.d05, 'master_cell_type', {4});
% d.d05.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, [d.d05s.offM]));
% clear cell_list_map;

d.d14 = load_data([piece '/data014'], loadopts);
d.d14 = read_stim_lisp_output(d.d14, ':2011-06-24-6:2011-06-24-6_f08_allcones');
d.d14.stimulus = parse_stim_rgbs(d.d14.stimulus);
 
cell_list_map = map_ei(d.d08s, d.d14, 'master_cell_type', {4}, 'electrode_threshold', 4);
d.d14.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, [d.d08s.offM]));
clear cell_list_map;

conerun = d.d08s;
datarun = d.d14;
conerun.rgcs = conerun.offM;
datarun.rgcs = datarun.offM_fd08s;
datarun.stimulus.triggers = datarun.triggers;

% % Just before, stability is okay, interestingly...
% if ~isfield(d, 'd13')
%     d.d13 = load_data([piece '/data013'], staopts);
% end
% cell_list_map = map_ei(d.d08s, d.d13, 'master_cell_type', {4});
% d.d13.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, [d.d08s.offM]));
% clear cell_list_map;

if ~isfield(d, 'd15')
    d.d15 = load_data([piece '/data015'], staopts);
end
cell_list_map = map_ei(d.d08s, d.d15, 'master_cell_type', {4});
d.d15.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, [d.d08s.offM]));
clear cell_list_map;

stablerun = [];
stablerun = d.d13;
stablerun.rgcs = d.d13.offM_fd08s;

allcones_plot_old(datarun, conerun, 'printpath', allconesprintpath, 'urgbs', datarun.stimulus.urgbs(1), 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});

% leave(keep_vars{:});


%% Allcones d14, try new allcones_plot
% Stability is terrible; totally unusable!

 
% d.d05 = load_data([piece '/data005'], loadopts);
% cell_list_map = map_ei(d.d05s, d.d05, 'master_cell_type', {4});
% d.d05.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, [d.d05s.offM]));
% clear cell_list_map;

d.d14 = load_data([piece '/data014'], loadopts);
d.d14 = read_stim_lisp_output(d.d14, ':2011-06-24-6:2011-06-24-6_f08_allcones');
d.d14.stimulus = convert_stimulus_to_combined_maps(d.d14.stimulus);
d.d14.stimulus = parse_stim_rgbs(d.d14.stimulus);
d.d14.stimulus.triggers = d.d14.triggers(1:2:end);

cell_list_map = map_ei(d.d08s, d.d14, 'master_cell_type', {4}, 'electrode_threshold', 4);
d.d14.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, [d.d08s.offM]));
clear cell_list_map;

conerun = d.d08s;
datarun = d.d14;
% conerun.rgcs = conerun.offM;
% datarun.rgcs = datarun.offM_fd08s;
datarun.stimulus.triggers = datarun.triggers;


if ~isfield(d, 'd15')
    d.d15 = load_data([piece '/data015'], staopts);
end
cell_list_map = map_ei(d.d08s, d.d15, 'master_cell_type', {4});
d.d15.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, [d.d08s.offM]));
clear cell_list_map;


stablerun = [];
stablerun = d.d15;
% stablerun.rgcs = d.d15.offM_fd08s;

conerun.rgcs = conerun.offM(1);
datarun.rgcs = datarun.offM_fd08s(1);
stablerun.rgcs = stablerun.offM_fd08s(1);
allcones_plot(conerun, datarun, datarun.stimulus.mapims{1}, 5, 'stablerun', stablerun, 'rfopts', {'scaled_up', 10});

% leave(keep_vars{:});


%% data005
d.d05 = load_data([piece '/data005'], staopts);

%% data009
d.d09 = load_data([piece '/data009'], loadopts);


%% data010
d10 = load_data('2011-06-24-6/data010');
d10 = load_neurons(d10);
d10 = load_ei(d10, 'all');


%% data011
d.d11 = load_data('2011-06-24-6/data011', loadopts);


%% data012
d12 = load_data('2011-06-24-6/data012');
d12 = load_neurons(d12);
d12 = load_ei(d12, 'all');


%% data013
d.d13 = load_data('2011-06-24-6/data013', staopts);


%% stimuli
s09 = read_stim_lisp_output('2011-06-24-6/stimuli/s09');
s10 = read_stim_lisp_output('2011-06-24-6/stimuli/s10');
s11 = read_stim_lisp_output('2011-06-24-6/stimuli/s11');
s12 = read_stim_lisp_output('2011-06-24-6/stimuli/s12');

s09.mapdir = [server_path '2011-06-24-6/stimuli/2011-06-24-6_data005_1and2'];
s10.mapdir = [server_path '2011-06-24-6/stimuli/2011-06-24-6_f05_12serial'];
s11.mapdir = [server_path '2011-06-24-6/stimuli/2011-06-24-6_f08_1and2'];
s12.mapdir = [server_path '2011-06-24-6/stimuli/2011-06-24-6_f08_12serial'];

mapstimuli = {s09 s10 s11 s12};
for s = 1:length(mapstimuli)
    stim = mapstimuli{s};
    mapims = {};
    mapnyc = {};
    for i = stim.mapindices
        fprintf('.');
        mapims{i+1} = dlmread([stim.mapdir sprintf('/map-%.4d.txt', i)]);
        mapnyc{i+1} = map2manhattan(mapims{i+1}, 'quiet', true);
    end; 
    mapstimuli{s}.mapims = mapims;
    mapstimuli{s}.mapnyc = mapnyc;
    fprintf('\n');
end; clear s i stim mapims;

s09 = mapstimuli{1};
s10 = mapstimuli{2};
s11 = mapstimuli{3};
s12 = mapstimuli{4};
clear mapstimuli;


%% Map d09 to d05s on EIs, get handpicked off midgets
cell_list_map = map_ei(d.d05s, d.d09, 'master_cell_type', {4});
d.d09.offMfd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
clear cell_list_map;

% The upshot
d.d09.offMfd05s = {50 1163 528 7206 2986 6676 4666 5656 5989};


%% Map d10 to d05str on EIs, get handpicked off midgets
cell_list_map = map_ei(d05str, d10, 'master_cell_type', {4});
off_midgets_fd05str_d10 = cell_list_map(get_cell_indices(d05str, off_midgets_d05str));
clear cell_list_map;

% The upshot
off_midgets_fd05str_d10 = {[] [] [] [] [] [] [] 5389 5916};


%% Map d11 to d08s on EIs, get handpicked off parasols
cell_list_map = map_ei(d.d08s, d.d11, 'master_cell_type', {2});
d.d11.offPfd08s_d11 = cell_list_map(get_cell_indices(d.d08s, d.d08s.offP));
clear cell_list_map;

% The upshot; manually added 7595
d.d11.offPfd08s = {5392 1381 1927 7596 6361 3303};


%% Map d12 to d08str on EIs, get handpicked off parasols
cell_list_map = map_ei(d08str, d12, 'master_cell_type', {2});
off_parasols_fd08str_d12 = cell_list_map(get_cell_indices(d08str, off_parasols_d08str));
clear cell_list_map;

% The upshot
off_parasols_fd08str_d12 = {5462 1396 1880 [] 6421 3286};


%% Searching d09


%% Build additivity structs
additivity_d05.conerun = d.d05s;
additivity_d05.conerun.addrgcs = d.d05s.offM;
additivity_d05.addrun  = d.d09;
additivity_d05.addstim = s09;
additivity_d05.addrun.addrgcs  = d.d09.offMfd05s;
additivity_d05.stablerun = d.d13;


additivity_d08.conerun = d08str;
additivity_d08.conerun.addrgcs = off_parasols_d08str;
additivity_d08.addrun  = d11;
additivity_d08.addstim = s11;
additivity_d08.addrun.addrgcs  = off_parasols_fd08str_d11;
additivity_d08.stablerun = d13;


%% Plot rasters d09: superposition, contrast-response
additivity_compound_plot(additivity_d05);
% These all look pretty good; the stability run seems like it has a few too
% many cones; might want to raise magic number and rerun?


%% Plot rasters d10: crosstalk control
crosstalk_raster_plot(d10, s10, off_midgets_fd05str_d10);


%% Plot rasters d11: superposition, contrast-response
additivity_compound_plot(additivity_d08, 'all');
% These look pretty good, see d05 comments above.


%% Plot rasters d12: crosstalk control

% THESE ARE WRONG!
% Unfortunately, we ran the wrong set of maps; in the stimuli.lisp file it
% is clear that we ran the maps from 2011-06-24-6_f05_12serial instead of
% 2011-06-24-6_f08_12serial.  :(

% To see the cones that were actually being stimulated, run this:
map = s10.mapims{1};
for i = 2:length(s12.mapims)
    map = map + i*s10.mapims{i};
end
d = d08str;
cells = off_parasols_d08str;
scale = 2;
ax1 = imagesc(map(1:scale:end, 1:scale:end)); colormap gray
hold on
plot_rf_summaries(d, cells, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b')
axis equal; axis tight;
lock;
clear map d cells scale ax1;

% Now you can see the funky results, which are reasonable given the above
crosstalk_raster_plot(d12, s12, off_parasols_fd08str_d12);


%% Plot rasters d12: crosstalk control (using as control for d09 now, since the stimuli got mixed up)

% Map d12 to d05str on EIs, get handpicked off parasols
cell_list_map = map_ei(d05str, d12, 'master_cell_type', {8});
off_midgets_fd05str_d12 = cell2mat(cell_list_map(get_cell_indices(d05str, off_midgets_d05str)));
clear cell_list_map;

% The upshot
off_midgets_fd05str_d12 = [52 1293 647 7202 2986 5988];

%crosstalk_raster_plot(d12, s12, off_midgets_fd05str_d12, 'maps', [0:8; 9:17]);


%% Contrast response analysis
cr09 = additivity_contrast_response({d09 d10 d12}, {s09 s10 s12}, {off_midgets_fd05str_d09 off_midgets_fd05str_d10 off_midgets_fd05str_d12}, 0.75, 0.02, 'start', -0.1);
cr11 = additivity_contrast_response({d11}, {s11}, {off_parasols_fd08str_d11}, 0.75, 0.02, 'start', -0.1);
cr09(:,4,:) = cr09(:,1,:) + cr09(:,2,:);
cr11(:,4,:) = cr11(:,1,:) + cr11(:,2,:);


%% Compare midget C/Rs scaled by STA weights to least squares
crstruct09.crfs = cr09;
crstruct09.crrun = d09;
crstruct09.crstim = s09;
crstruct09.crcells = off_midgets_fd05str_d09;
crstruct09.conerun = d05str;
crstruct09.conecells = off_midgets_d05str;
crstruct09.stablerun = d13; % Should be able to change this to d08
inds = {[1 2] [3 4] [5 6] [7 8] [9]};
for i = 1:length(inds)
    figure;
    cr_compound_plot(crstruct09, inds{i});
end; clear i inds


%% Compare parasol C/Rs scaled by STA weights to least squares
crstruct11.crfs = cr11;
crstruct11.crrun = d11;
crstruct11.crstim = s11;
crstruct11.crcells = off_parasols_fd08str_d11;
crstruct11.conerun = d08str;
crstruct11.conecells = off_parasols_d08str;
crstruct11.stablerun = d13;
inds = {[1 2] [3 4] [5 6]};
for i = 1:length(inds)
    figure;
    cr_compound_plot(crstruct11, inds{i});
end; clear i inds


%% Compare sensitivity
subplot(9, 1, 1:5);
hm = bar(0:25:250, [histc(squeeze(cr09(1,1,:)), 0:25:250) histc(squeeze(cr09(1,2,:)), 0:25:250)], 0.95, 'stacked'); % Midgets
set(gca, 'XTick', [], 'YTick', 0:5);
axis([0 250 0 5.25]);
legend(hm, {'1st' '2nd'});

subplot(9, 1, 7:9);
hp = bar(0:25:250, [histc(squeeze(cr11(1,1,:)), 0:25:250) histc(squeeze(cr11(1,2,:)), 0:25:250)], 0.95, 'stacked'); % Parasols
set(gca, 'XTick', 0:50:250, 'YTick', 0:3);
axis([0 250 0 3.25]);

clear hp hm;


%% Plot map with STA fits over
map = s11.mapims{3};
d = d08str;
cells = off_parasols_d08str;

scale = 2;
ax1 = imagesc(map(1:scale:end, 1:scale:end)); colormap gray
hold on

plot_rf_summaries(d, cells, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b')
axis equal; axis tight;

clear map d cells scale ax1;


%% Plot STA with voronoi over
% d = d05str;
% cones = d05str.cones;
% cell = off_midgets_d05str(2);
% % plot_sta(d, cell);
% plot_rf(d, cell, 'fit', false);
% 
% hold on;
% [vx, vy] = voronoi(cones.centers(:,1), cones.centers(:,2));
% h = plot(vx, vy, 'b');
% set(h(1:end-1), 'xliminclude', 'off', 'yliminclude', 'off')
% 
% autozoom_to_fit(d, cell);
% 
% clear d cones cell vx vy h;


%% Plot map with voronoi over
% map = s09.mapims{3};
% cones = d05str.cones;
% voronoi_masks_scale = 2;
% imagesc(map); hold on; voronoi(cones.centers(:,1) .* voronoi_masks_scale, cones.centers(:,2) .* voronoi_masks_scale);
% axis equal; axis tight;
% 
% clear map cones voronoi_masks_scale;


%% Good RF On M
d.d03 = load_data([piece '/data003'], staopts);
d.d04 = load_data([piece '/data004'], staopts);
d.d05 = load_data([piece '/data005'], staopts);
d.d03.mapd04 = map_ei(d.d04, d.d03);
d.d05.mapd04 = map_ei(d.d04, d.d05);

d.d04.rgc = 1863;
d.d03.rgc = d.d03.mapd04(get_cell_indices(d.d04, d.d04.rgc));
d.d05.rgc = d.d05.mapd04(get_cell_indices(d.d04, d.d04.rgc));

d.d03.rf = get_rf(d.d03, d.d03.rgc{1});
d.d04.rf = get_rf(d.d04, d.d04.rgc);
d.d05.rf = get_rf(d.d05, d.d05.rgc{1});
rfcompound = d.d03.rf + sum(d.d04.rf, 3) + sum(d.d05.rf, 3);
datarun = d.d04;
datarun.stas.rfs{get_cell_indices(datarun, d.d04.rgc)} = rfcompound;
plot_rf(datarun, d.d04.rgc, 'fit', false, 'autozoom', true, 'scale', 10);


%% Good RF Off M
d.d04 = load_data([piece '/data004'], staopts);
d.d05 = load_data([piece '/data005'], staopts);
d.d04.mapd05 = map_ei(d.d05, d.d04);

d.d05.rgc = 51;
d.d04.rgc = d.d04.mapd05(get_cell_indices(d.d05, d.d05.rgc));

d.d04.rf = get_rf(d.d04, d.d04.rgc{1});
d.d05.rf = get_rf(d.d05, d.d05.rgc);
rfcompound = d.d04.rf + d.d05.rf;
datarun = d.d05;
datarun.stas.rfs{get_cell_indices(datarun, d.d05.rgc)} = rfcompound;
plot_rf(datarun, d.d05.rgc, 'fit', false, 'autozoom', true, 'scale', 10, 'color_transform', ones(3));


%% Good RF Off P
% Not as good as 2012-09-13-2
%
% d.d04 = load_data([piece '/data004'], staopts);
% d.d05 = load_data([piece '/data005'], staopts);
% d.d04.mapd05 = map_ei(d.d05, d.d04);
% 
% d.d04.rgc = 5881;
% d.d05.rgc = 5851;
% 
% d.d04.rf = get_rf(d.d04, d.d04.rgc);
% d.d05.rf = get_rf(d.d05, d.d05.rgc);
% rfcompound = d.d04.rf + d.d05.rf;
% datarun = d.d05;
% datarun.stas.rfs{get_cell_indices(datarun, d.d05.rgc)} = rfcompound;
% plot_rf(datarun, d.d05.rgc, 'fit', false, 'autozoom', true, 'scale', 10, 'color_transform', ones(3));


%% Cone projective field - d04/d09

% No ON Parasol matches
d.d09.mapd04 = map_ei(d.d04, d.d09);
additivity_raster_plot(d.d09, s09, d.d09.mapd04(get_cell_indices(d.d04, {1})), 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);

% Off parasol
cell_list_map = map_ei(d.d04, d.d09, 'master_cell_type', {2});
d.d09.offP = cell_list_map(get_cell_indices(d.d04, {2}));
clear cell_list_map;
additivity_raster_plot(d.d09, s09, d.d09.offP, 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);


%% Cone projective field - d05s/d09

% First lets take a look through all the off midgets as a control
cell_list_map = map_ei(d.d05s, d.d09, 'master_cell_type', {4});
all_off_midgets_fd05s_d09 = cell_list_map(get_cell_indices(d.d05s, {4}));
clear cell_list_map;
additivity_raster_plot(d.d09, s09, all_off_midgets_fd05s_d09, 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);

% Only three that weren't handpicked: 1732, 5240, 5657
% Only two with real responses, but these are pretty big!
%   1732 (d05str: 96), 5240 (d05str: 106)
% Interesting; d05str_cell96 looks weird.  d05str_cell106 looks like real
% cone overlap; is it really a midget?



% Off parasol
cell_list_map = map_ei(d.d05s, d.d09, 'master_cell_type', {2});
all_off_parasols_fd05s_d09 = cell_list_map(get_cell_indices(d.d05s, {2}));
clear cell_list_map;
additivity_raster_plot(d.d09, s09, all_off_parasols_fd05s_d09, 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);

% Very nice ones:
%           242, 1247, 1732, 1819, 4430, 4598, 5043, 5736, 5777, 6391, 6422, 7340
% d05str:    49,   75,   98,  150,   16,   23,  103,  107,   30,  111,  115,  147


% No on midgets :(



% On parasol; no matches :(
cell_list_map = map_ei(d.d05s, d.d09, 'master_cell_type', {1});
all_on_parasols_fd05s_d09 = cell_list_map(get_cell_indices(d.d05s, {1}));
clear cell_list_map;
additivity_raster_plot(d.d09, s09, all_on_parasols_fd05s_d09, 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);


%% Cone projective field - d05/d09
d.d09.mapd05 = map_ei(d.d05, d.d09);
additivity_raster_plot(d.d09, s09, d.d09.mapd05(get_cell_indices(d.d05, {1})), 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);
% No ON Parasol matches

%% Cone projective field d05str/d09; plotting
plot_rf_summaries(d05str, off_midgets_d05str, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k')
plot_rf_summaries(d05str, [49 75 98 150 16 23 103 107 30 111 115 147], 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b')

% Okay, but parasol is not so pretty
additivity_raster_plot(d09, s09, [242 50], 'intensities', [-0.44], 'figure_by_cell_id', true); % Better
cone2 = get_cones_by_weight_rank(d05str, 50, [2]);
plot_voronoi_over_rf(d05str, 50, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0]); lock
%plot_voronoi_over_rf(d05str, 49, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0]); lock
d05.cones = d05str.cones; plot_voronoi_over_rf(d05, 422, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0]); lock

% New plotting
plot_rf(d05str, 50, 'fit', false, 'scale', 10);
scalemap = 1 / d05str.cones.mosaic.voronoi_masks_scale;
patchmanhattan(s09.mapnyc{2}, 'scale', [scalemap scalemap], 'colors', [0 1 0], 'patchopts', {'LineWidth', 2});
autozoom_to_fit(d05str, 50, 6, 1, 1)

plot_rf(d05, 422, 'fit', false, 'scale', 10);
patchmanhattan(s09.mapnyc{2}, 'scale', [scalemap scalemap], 'colors', [0 1 0], 'patchopts', {'LineWidth', 2});
autozoom_to_fit(d05, 422, 4, 1, 1)
clear scalemap;

% % Same problem
% additivity_raster_plot(d09, s09, [4598 4666], 'intensities', [-0.44 0.44], 'figure_by_cell_id', true); % Better
% cone1 = get_cones_by_weight_rank(d05str, 15, [1]);
% plot_voronoi_over_rf(d05str, 15, 'highlight_cones', cone1, 'highlight_rgba', [0 1 0]); lock
% plot_voronoi_over_rf(d05str, 23, 'highlight_cones', cone1, 'highlight_rgba', [0 1 0]); lock
% 
% % Cones pretty far from parasol...
% additivity_raster_plot(d09, s09, [6422 5989], 'intensities', [-0.44 0.44], 'figure_by_cell_id', true); % Better
% cones12 = get_cones_by_weight_rank(d05str, 42, [1 2]);
% plot_voronoi_over_rf(d05str, 42, 'highlight_cones', cones12, 'highlight_rgba', [1 0 0; 0 1 0]); lock
% plot_voronoi_over_rf(d05str, 115, 'highlight_cones', cones12, 'highlight_rgba', [1 0 0; 0 1 0]); lock


% Usable
additivity_raster_plot(d09, s09, [7340 7206], 'intensities', [-0.44 0.44], 'figure_by_cell_id', true); % Better
cone2 = get_cones_by_weight_rank(d05str, 137, [2]);
plot_voronoi_over_rf(d05str, 137, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0]); lock
%plot_voronoi_over_rf(d05str, 147, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0]); lock
% Hack this to use better STA from d05 post analysis!!
d05.cones = d05str.cones; plot_voronoi_over_rf(d05, 316, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0]); lock


%% Cone projective field - d08s/d11

% First lets take a look through all the off parasols as a control
cell_list_map = map_ei(d08str, d11, 'master_cell_type', {2});
all_off_parasols_fd08s_d11 = cell_list_map(get_cell_indices(d05str, {2}));
clear cell_list_map;
additivity_raster_plot(d11, s11, all_off_parasols_fd08s_d11, 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);



% Off midget
cell_list_map = map_ei(d08str, d11, 'master_cell_type', {4});
all_off_midgets_fd08s_d11 = cell_list_map(get_cell_indices(d05str, {4}));
clear cell_list_map;
additivity_raster_plot(d11, s11, all_off_midgets_fd08s_d11, 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);
% Not good



% On parasol
cell_list_map = map_ei(d08str, d11, 'master_cell_type', {1});
all_on_parasols_fd08s_d11 = cell_list_map(get_cell_indices(d05str, {1}));
clear cell_list_map;
additivity_raster_plot(d11, s11, all_on_parasols_fd08s_d11, 'intensities', [-0.44 0.44], 'figure_by_cell_id', true);
% 1382 (d08str_cell46) interesting, but wrong direction...  Pretty close to
% edge.


%% Find specific matches from mapping
mapped = all_off_parasols_fd05s_d09;
from = d05str.cell_types{2}.cell_ids;
interest = [242, 1247, 1732, 1819, 4430, 4598, 5043, 5736, 5777, 6391, 6422, 7340];

for i = 1:length(mapped)
    cell = mapped{i};
    if ~isempty(intersect(interest, cell))
        disp([num2str(cell) ': ' num2str(from(i))]);
    end
end

clear i mapped from interest cell ans


%% Compare additivity DEPRECATED
% figure; 
% subplot(1, 2, 1); 
% hold on;
% plot(-10:400, -10:400, 'k--');
% h201106246_d09 = plot(squeeze(cr09(:,4,:)), squeeze(cr09(:,3,:)), 'b.');
% set(h201106246_d09, 'MarkerSize', 20);
% axis equal tight;
% title 'Off Midgets, 2011-06-24-6';
% ylabel 'Observed Firing Rate (Hz)';
% xlabel 'Linear Prediction (Hz)';
% set(gca, 'XTick', 0:100:400);
% set(gca, 'YTick', 0:100:400);
% 
% subplot(1, 2, 2);
% hold on;
% plot(-10:400, -10:400, 'k--');
% h201106246_d11 = plot(squeeze(cr11(:,4,:)), squeeze(cr11(:,3,:)), 'b.');
% set(h201106246_d11, 'MarkerSize', 20);
% axis equal tight;
% title 'OFF Parasols, 2011-06-24-6'
% xlabel 'Linear Prediction (Hz)';
% set(gca, 'YTick', []);
% set(gca, 'XTick', 0:100:400);
% 
% linkaxes(get(gcf, 'Children'));