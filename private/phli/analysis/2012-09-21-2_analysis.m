% Basics
piece = '2012-09-21-2';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%% Streamed runs and cell picks
d.d09s = load_data([piece '/streamed/data009/data009'], staopts);
d.d09s.offM = [151 181 541 691 1051 1276 1351 2658 2672 2731 4141 4216 4381 5086 5491 6046 6346 6406 7006 7696];

d.d11s = load_data([piece '/streamed/data011/data011'], staopts);
d.d11s.offM = [152 181 542 691 1051 1276 1352 2673 2672 2731 4141 3481 4471 5086 5491 6046 6346 6406 7006 7696];
%d.d11s.onPsearched = [2116 4456]; Didn't really pan out
d.d11s.onPsurround = [5568];


%% data013 allcones09 analysis
% Got all but one mapped, a few more missing stability.  Can check
% concatenated for the rest if worthwhile.
allconesprintpath = printpath('allcones', piece);

if ~isfield(d, 'd13'), d.d13 = load_data([piece '/data013'], loadopts); end
d.d13 = read_stim_lisp_output(d.d13);

cell_list_map = map_ei(d.d09s, d.d13, 'master_cell_type', 'all');
d.d13.offM_f09s = cell_list_map(get_cell_indices(d.d09s, d.d09s.offM));
clear cell_list_map;

conerun = d.d09s;
datarun = d.d13;
conerun.rgcs = d.d09s.offM;
datarun.rgcs = d.d13.offM_f09s;
datarun.stimulus.triggers = datarun.triggers(1:2:end);

if ~isfield(d, 'd14'), d.d14 = load_data([piece '/data014'], staopts); end
cell_list_map = map_ei(d.d09s, d.d14, 'master_cell_type', 'all');
d.d14.offM_f09s = cell_list_map(get_cell_indices(d.d09s, d.d09s.offM));
clear cell_list_map;

stablerun = d.d14;
stablerun.rgcs = d.d14.offM_f09s;

allcones_plot(datarun, conerun, 'coneind', 'localmax', 'urgbs', {[1 1 1].*-0.336 [1 1 1].*-0.48}, 'printpath', allconesprintpath, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});

leave(keep_vars);


%% data015cm C/R
crprintpath = printpath('cr', piece);

if ~isfield(d, 'd15cm'), d.d15cm = load_data([piece '/d14-15-16-norefit/data015/data015'], loadopts); end
if ~isfield(d, 'd16cm'), d.d16cm = load_data([piece '/d14-15-16-norefit/data016/data016'], staopts); end

% Map from d11s to concatenated (valid also for individuals from concatenated since norefit)
if ~isfield(d.d15cm, 'mapd11s')
    if ~isfield(d, 'd14_15_16_n'), d.d14_15_16_n = load_data([piece '/d14-15-16-norefit/d14-15-16-norefit'], loadopts); end
    d.d15cm.mapd09s = map_ei(d.d09s, d.d14_15_16_n);
    d.d15cm.mapd11s = map_ei(d.d11s, d.d14_15_16_n);
end

conerun = d.d11s;
conerun.rgcs = [conerun.offM];
rasterrun = d.d15cm;
rasterrun.rgcs = rasterrun.mapd11s(get_cell_indices(conerun, conerun.rgcs));
stablerun = d.d16cm;
stablerun.rgcs = rasterrun.rgcs;
stablerun = set_polarities(stablerun, 'cell_specs', {{8}}, 'polarities', -1);

[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d15cm d.d11s d.d16cm] = deal(cleared{:});
% Some pretty nice ones here, but the movement is bad.  Some of the
% IDs seem to be wrong; no regions over the cell?
%   See note on U/D below


%% data015cm U/D
udprintpath = printpath('ud', piece);

if ~isfield(d, 'd15cm'), d.d15cm = load_data([piece '/d14-15-16-norefit/data015/data015'], loadopts); end
if ~isfield(d, 'd16cm'), d.d16cm = load_data([piece '/d14-15-16-norefit/data016/data016'], staopts); end

% Map from d11s to concatenated (valid also for individuals from concatenated since norefit)
if ~isfield(d.d15cm, 'mapd11s')
    if ~isfield(d, 'd14_15_16_n'), d.d14_15_16_n = load_data([piece '/d14-15-16-norefit/d14-15-16-norefit'], loadopts); end
    d.d15cm.mapd11s = map_ei(d.d11s, d.d14_15_16_n);
end
d.d16cm.mapd11s = d.d15cm.mapd11s;

d.d15cm = read_stim_lisp_output(d.d15cm);
d.d15cm.stimulus = parse_stim_rgbs(d.d15cm.stimulus);

conerun = d.d11s;
conerun.wwrgcs = conerun.offM;

stablerun = d.d16cm;
stablerun.wwrgcs = stablerun.mapd11s(get_cell_indices(conerun, conerun.wwrgcs));
stablerun = set_polarities(stablerun, 'cell_specs', {{1} {8}});

ud15.wwrun = d.d15cm;
ud15.wwrun.wwrgcs = ud15.wwrun.mapd11s(get_cell_indices(conerun, conerun.wwrgcs));
ud15.conerun = conerun;
ud15.stablerun = stablerun;
ww_compound_plot(ud15, 'all', 'triggers', ud15.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear conerun stablerun udprintpath

% No stimulation and in a different cell type...
% 181, 2731, 4141, 5086, 6046


% For J Freeman subunits analysis
save_udraster_data(ud15, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));



%% data015cm searching
if ~isfield(d, 'd15cm'), d.d15cm = load_data([piece '/d14-15-16-norefit/data015/data015'], loadopts); end
if ~isfield(d, 'd14_15_16_n'), d.d14_15_16_n = load_data([piece '/d14-15-16-norefit/d14-15-16-norefit'], loadopts); end
map_d15cm_d11s = map_ei(d.d11s, d.d14_15_16_n);
d.d15cm = read_stim_lisp_output(d.d15cm);
d.d15cm.stimulus = parse_stim_rgbs(d.d15cm.stimulus);
urgb = d.d15cm.stimulus.urgb;

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d11s;
rasterrun = d.d15cm;
map = map_d15cm_d11s;

% Search all ON parasols
celltype = 1;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Number of interesting ones, vague hope for single cone on 2116, 4456 (didn't really pan out)
% Single cone surround possibility: 2116, 5568 (2116 too messy; can still test other)

%% Stability check, Jeremy Freeman's repeat analysis
d.d07 = load_data([piece '/data007'], staopts);
d.d07 = load_cones(d.d07);
d.d09 = load_data([piece '/data009'], staopts);
d.d09 = load_cones(d.d09);
overlay_cone_mosaics(d.d07, d.d09);
plot_rf_summaries(d.d07, {4}, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y');