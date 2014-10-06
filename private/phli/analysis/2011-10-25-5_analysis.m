%% Basic load data
piece = '2011-10-25-5';
udprintpath = [printpath() '/ud'];

d01s = load_data([piece '/streamed/data001-0/data001-0']);
d04s = load_data([piece '/streamed/data004-0/data004-0']);
d03 = load_data([piece '/data003']);
d04 = load_data([piece '/data004']);
d05 = load_data([piece '/data005']);
d06 = load_data([piece '/data006']);
d07 = load_data([piece '/data007']);

dataruns = {'d01s' 'd04s' 'd03' 'd04' 'd05' 'd06' 'd07'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% EIs
dataruns = {'d01s' 'd04s' 'd03' 'd04' 'd05' 'd06' 'd07'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_ei(datarun, 'all');
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% STAs and cones
dataruns = {'d01s' 'd04s' 'd04' 'd07'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = set_polarities(datarun);
    datarun = get_sta_fits_from_vision(datarun);
    
    datarun = load_cones(datarun, 1);
    datarun = make_mosaic_struct(datarun);
    datarun.cones.mosaic.voronoi_masks_scale = 2;
    datarun = make_voronoi_masks(datarun);
    
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i

% Fix these
d04s = load_cones(d04s, 'acquisition');
d07.cones.mosaic.voronoi_masks_scale = 3;
d07 = make_voronoi_masks(d07);


%% Picks
% 2011-10-25-5 om fd01
offmidgets_d01s = [68 184 698 1491 2554 2583 3169 6408 6631 454 948 5222];
offparasols_d04s = [1985 1277 512 6046 7171 2176];


%% Stimuli
d03 = read_stim_lisp_output(d03);
d03.stimulus = parse_stim_rgbs(d03.stimulus);
d06 = read_stim_lisp_output(d06);
d06.stimulus = parse_stim_rgbs(d06.stimulus);

%% Map d03 to d01s on EIs, get handpicked off midgets
cell_list_map = map_ei(d01s, d03, 'master_cell_type', {4});
offmidgets_fd01s_d03 = cell_list_map(get_cell_indices(d01s, offmidgets_d01s));
clear cell_list_map;


%% Map d04s to d01s on EIs, get handpicked off midgets
cell_list_map = map_ei(d01s, d04s, 'master_cell_type', {4});
offmidgets_fd01s_d04s = cell_list_map(get_cell_indices(d01s, offmidgets_d01s));
clear cell_list_map;


%% Map d05 to d01s on EIs, get handpicked off midgets
cell_list_map = map_ei(d01s, d05, 'master_cell_type', {4});
offmidgets_fd01s_d05 = cell_list_map(get_cell_indices(d01s, offmidgets_d01s));
clear cell_list_map;


%% Map d06 to d04s on EIs, get handpicked off parasols
cell_list_map = map_ei(d04s, d06, 'master_cell_type', {2});
offparasols_fd04s_d06 = cell_list_map(get_cell_indices(d04s, offparasols_d04s));
clear cell_list_map;


%% Map d07 to d04s on EIs, get handpicked off parasols
cell_list_map = map_ei(d04s, d07, 'master_cell_type', {2});
offparasols_fd04s_d07 = cell_list_map(get_cell_indices(d04s, offparasols_d04s));
clear cell_list_map;


%% Searching d03
conerun = d01s;
rasterrun = d03;
map = map_ei(d01s, d03);

urgb = rasterrun.stimulus.urgb;
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};

% Not much for ON Parasols

% Search all OFF Parasols
celltype = 2;
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% 6151, 6332 pretty good

% Search all M+
celltype = 3;
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% 2103 nice RF and has surround response



%% U/D d03, new for subunits paper

conerun = d01s;
conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
conerun = set_polarities(conerun);
conerun = get_sta_fits_from_vision(conerun);
conerun = load_cones(conerun, 'acquisition');
stablerun = d04s;
stablerun = load_sta(stablerun, struct('load_sta', [], 'guess_stimulus', false));
stablerun = set_polarities(stablerun);
stablerun = get_sta_fits_from_vision(stablerun);
rasterrun = d03;

figure
outstructs = ud_compound_plot(conerun, rasterrun, 6408, 'triggers', rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rfopts', {'scaled_up', 10}, 'rasteropts', {'axopts', {{'Box', 'off', 'XTick', [0 0.6]}}});
outstructs = cell2mat(invertselect(outstructs, @isempty));
AX = vertcat(outstructs.AX);
set(AX(:,1), 'YTick', []);
set(AX(:,2), 'YTick', [0 40]);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(:,1), 'XTickLabel', []);
set(AX(setdiff(1:20, 4:4:20),2), 'XTickLabel', []);
set(gcf, 'Position', [2148 69 1185 846]);

figure
outstructs = ud_compound_plot(conerun, rasterrun, 2583, 'triggers', rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rfopts', {'scaled_up', 10}, 'rasteropts', {'axopts', {{'Box', 'off', 'XTick', [0 0.6]}}});
outstructs = cell2mat(invertselect(outstructs, @isempty));
AX = vertcat(outstructs.AX);
set(AX(:,1), 'YTick', []);
set(AX(:,2), 'YTick', [0 40]);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(:,1), 'XTickLabel', []);
set(AX(setdiff(1:20, 4:4:20),2), 'XTickLabel', []);
set(gcf, 'Position', [2148 69 1185 846]);

figure
outstructs = ud_compound_plot(conerun, rasterrun, 948, 'triggers', rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rfopts', {'scaled_up', 10}, 'rasteropts', {'axopts', {{'Box', 'off', 'XTick', [0 0.6]}}});
outstructs = cell2mat(invertselect(outstructs, @isempty));
AX = vertcat(outstructs.AX);
set(AX(:,1), 'YTick', []);
set(AX(:,2), 'YTick', [0 40]);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(:,1), 'XTickLabel', []);
set(AX(setdiff(1:20, 4:4:20),2), 'XTickLabel', []);
set(gcf, 'Position', [2148 69 1185 846]);


%%
ww03.wwrun = d03;
ww03.wwrun.wwrgcs = offmidgets_fd01s_d03;
ww03.conerun = d01s;
ww03.conerun.wwrgcs = offmidgets_d01s;
ww03.stablerun = d04s;
ww03.stablerun.wwrgcs = offmidgets_fd01s_d04s;
ww_compound_plot(ww03, 'all', 'triggers', d03.triggers(1:2:end), 'printpath', udprintpath);

% For J Freeman subunits analysis
save_udraster_data(ww03, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


%% Separate plotting for figures
ww_compound_plot(ww03, 10, 'triggers', d03.triggers(1:2:end), 'scaled_up', 10, 'MarkerSize', 15, 'az_pad_factor', 4.5);
ww_compound_plot(ww03, 11, 'triggers', d03.triggers(1:2:end), 'scaled_up', 10, 'MarkerSize', 15, 'az_pad_factor', 4.5);

% figure;
% plot_voronoi_over_rf(ww03.conerun, ww03.conerun.wwrgcs(8), 'az_aspect_ratio', 1, 'scaled_up', 10);
% scalemap = 1 / ww03.conerun.cones.mosaic.voronoi_masks_scale;
% patchmanhattan(ww03.wwstim.mapnyc{1}, 'scale', [scalemap scalemap], 'colors', [1 0 0; 0 0.5 0; 0 0 1; 1 0 1], 'patchopts', {'LineWidth', 2});
% clear scalemap;

% 3 is a little interesting for cones 3/4; asymmetric cones 1/2 sketchy; matches instability
% 9 looks a little interesting


%%
ww06.wwrun = d06;
ww06.wwrun.wwrgcs = offparasols_fd04s_d06;
ww06.conerun = d04s;
ww06.conerun.wwrgcs = offparasols_d04s;
ww06.stablerun = d07;
ww06.stablerun.wwrgcs = offparasols_fd04s_d07;
ww_compound_plot(ww06, 'all', 'triggers', d06.triggers(1:2:end), 'printpath', udprintpath);

% fig 4, cell 4, id 6046 is quite pretty


% For J Freeman subunits analysis
save_udraster_data(ww06, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));



%% Separate plotting for figures

ww_compound_plot(ww06, 4, 'triggers', d06.triggers(1:2:end), 'scaled_up', 10, 'MarkerSize', 15);
ww_compound_plot(ww06, 5, 'triggers', d06.triggers(1:2:end), 'scaled_up', 10, 'MarkerSize', 15);

% figure;
% plot_rf(ww06.conerun, ww06.conerun.wwrgcs(4), 'scale', 10, 'fit', false);
% scalemap = 1 /ww06.conerun.cones.mosaic.voronoi_masks_scale;
% patchmanhattan(ww06.wwstim.mapnyc{1}, 'scale', [scalemap scalemap], 'colors', [1 0 0; 0 0.5 0; 0 0 1; 1 0 1], 'patchopts', {'LineWidth', 2});
% autozoom_to_fit(ww06.conerun, ww06.conerun.wwrgcs(4), 4, 1, 1);
% clear scalemap


%% Testing cone refinding

% Map d04 to d01s on EIs
cell_list_map = map_ei(d01s, d04, 'master_cell_type', {4});
offmidgets_fd01s_d04 = cell_list_map(get_cell_indices(d01s, offmidgets_d01s));
clear cell_list_map;

% Plot an RF and the overlaid cone runs with cone labels
plot_rf(d01s, offmidgets_d01s(2))
overlay_cone_mosaics(d01s, d04, 'clear', false, 'label', true)