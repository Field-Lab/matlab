%% Basics
piece = '2012-08-21-0';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'udprintpath'; 'allconesprintpath'; 'keep_vars'; 'd'};


%% Streamed runs and cell picks

d.d01s = load_data([piece '/streamed/data001/data001'], staopts);
d.d01s.offM = [31 513 693 1009 1788 1996 2973 3153 4399 4703 5431 6724 6953 1101 4338 4027];

d.d01s.offMcones = [504, 561, 533, 547
    411, 450, 456, 480
    176, 203, 195, 221
    530, 531, 542, 557
    292, 343, 371, 400
    323, 347, 369, 386
    865, 878, 887, 858
    1021, 1061, 1081, 1095
    1214, 1130, 1120, 1169
    1457, 1464, 1513, 1545
    1791, 1821, 1820, 1850
    1017, 1043, 1085, 1131
    786, 820, 821, 822
    67, 71, 81, 101
    1711, 1748, 1755, 1767
    1435, 1473, 1482, 1526];


%% data005 C/R
crprintpath = [];%printpath('cr', piece);

if ~isfield(d, 'd05'), d.d05 = load_data([piece '/data005'], loadopts); end
if ~isfield(d, 'd06'), d.d06 = load_data([piece '/data006'], staopts); end
if ~isfield(d.d05, 'offM_fd01s')
    cell_list_map = map_ei(d.d01s, d.d05, 'master_cell_type', {4});
    d.d05.offM_fd01s = cell_list_map(get_cell_indices(d.d01s, d.d01s.offM));
    clear cell_list_map;
end
if ~isfield(d.d06, 'offM_fd01s')
    cell_list_map = map_ei(d.d01s, d.d06, 'master_cell_type', {4});
    d.d06.offM_fd01s = cell_list_map(get_cell_indices(d.d01s, d.d01s.offM));
    clear cell_list_map;
end

conerun = d.d01s;
conerun.rgcs = conerun.offM;
rasterrun = d.d05;
rasterrun.rgcs = rasterrun.offM_fd01s;
stablerun = d.d06;
stablerun.rgcs = stablerun.offM_fd01s;

[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d05 d.d01s d.d06] = deal(cleared{:});



%% data005cm searching
if ~isfield(d, 'd01cm'), d.d01cm = load_data([piece '/d01_05_06_n/data001/data001'], staopts); end
map_d05cm_d01cm = num2cell(d.d01cm.cell_ids);

if ~isfield(d, 'd05cm'), d.d05cm = load_data([piece '/d01_05_06_n/data005/data005'], loadopts); end
d.d05cm = read_stim_lisp_output(d.d05cm);
d.d05cm.stimulus = parse_stim_rgbs(d.d05cm.stimulus);
urgb = d.d05cm.stimulus.urgb;

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d01cm;
rasterrun = d.d05cm;
map = map_d05cm_d01cm;

% Nothing for ON Parasols


%% data005 searching
if ~isfield(d, 'd05'), d.d05 = load_data([piece '/data005'], loadopts); end
map_d05_d01s = map_ei(d.d01s, d.d05);
d.d05 = read_stim_lisp_output(d.d05);
d.d05.stimulus = parse_stim_rgbs(d.d05.stimulus);
urgb = d.d05.stimulus.urgb;

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d01s;
rasterrun = d.d05;
map = map_d05_d01s;

% Nothing for ON parasols
% Nothing much for ON Midgets

% Search all OFF parasols
celltype = 2;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% A few pretty nice centers in here


%% Crest analysis
% d.d03 = load_data([piece '/data003']);


%% Null stimulus data008
% Copied from 2012-09-13-2; still working out the kinks.

if ~isfield(d, 'd06s'), d.d06s = load_data([piece '/streamed/data006/data006'], staopts); end
if ~isfield(d, 'd08'), d.d08 = load_data([piece '/data008'], loadopts); end
if ~isfield(d, 'd09s'), d.d09s = load_data([piece '/streamed/data009/data009'], staopts); end

% Have to be careful not to try to parse polylines because the 2nd stim region
% isn't well formed and this will crash.
d.d08 = read_stim_lisp_output(d.d08, [], false, false);
d.d08.stimulus = parse_stim_rgbs(d.d08.stimulus);

d.d08.mapd06s = map_ei(d.d06s, d.d08);
d.d09s.mapd06s = map_ei(d.d06s, d.d09s);

ploturgbs = {3 2 1; 6 5 4};
titles = {-0.48 -0.24 -0.12};
conerun = d.d06s;
rasterrun = d.d08;
map = d.d08.mapd06s;
stablerun = d.d09s;
stablemap = d.d09s.mapd06s;

celltype = 2;
%  5851 best (d06s)

for cellid = conerun.cell_types{celltype}.cell_ids(1:end)
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'add_width', 1, 'titles', titles);

    sanesubplot(2, 7, {1, 7});
    ax = plot_rf_stimmap(conerun, cellid, flipud(rasterrun.stimulus.mapnyc{1}), 'fit', false, 'az_pad_factor', 4, 'colors', [0.5 0.5 1; 0 0 1]);

    stableid = stablemap{cellnum};
    if ~isempty(stableid)
        sanesubplot(2, 7, {2, 7});
        ax2 = plot_rf_stimmap(stablerun, stableid, flipud(rasterrun.stimulus.mapnyc{1}), 'fit', false, 'autozoom', false, 'colors', [0.5 0.5 1; 0 0 1]);
        set(ax2, 'XLim', get(ax, 'XLim'), 'YLim', get(ax, 'YLim'));
    end
end

%% Stability check, Jeremy Freeman's repeat analysis
d.d01 = load_data([piece '/data001'], staopts);
d.d01 = load_cones(d.d01, 1);
d.d01s = load_data([piece '/streamed/data001/data001'], staopts);
d.d01s = load_cones(d.d01s, 2);
d.d06s = load_data([piece '/streamed/data006/data006'], staopts);
d.d06s = load_cones(d.d06s, 3);
overlay_cone_mosaics(d.d01s, d.d06s);
plot_rf_summaries(d.d01s, {4}, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y');