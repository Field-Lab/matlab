% Basics
piece = '2012-09-13-2';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};

% This one has some interesting large cells that have raster responses; see search at bottom

%% Streamed runs and cell picks
% Too much movement from 1 to 5 to 9 for combined STAs to look better than
% individual in most cases

% Movement, but some good ones
d.d01s = load_data([piece '/streamed/data001/data001'], staopts);
d.d01s = load_txt_cell_types(d.d01s, 'classification');
d.d01s.onM = [1909 2225 2267 2449 2795 4056 4337 5556 6079 6153 6409 6679 6802 6917 7292];

% Movement, but some good ones
d.d05s = load_data([piece '/streamed/data005/data005'], staopts);
d.d05s.offM = [1396 1607 1727 1816 1966 2026 2116 2386 2942 6991];
d.d05s.offP = [46 376 1231 4441 5371 5551 6496];

% Moved too much probably
d.d09s = load_data([piece '/streamed/data009/data009'], staopts);
d.d09s.offM = [151 1396 1741 1816 1863 1966 2026 2116 6796 6991];
d.d09s.offP = [376 1471 4636 4771 5551];


%% Good RF off P
d.d05 = load_data([piece '/data005'], staopts);
plot_rf(d.d05, 6031, 'fit', false, 'autozoom', true, 'scale', 10);


%% data004 U/D onM
udprintpath = [];

d.d01 = load_data([piece '/data001'], loadopts);
d.d04 = load_data([piece '/data004'], loadopts);
d.d05 = load_data([piece '/data005'], loadopts);

cell_list_map = map_ei(d.d01s, d.d01, 'master_cell_type', {3});
d.d01.onM_fd01s = cell_list_map(get_cell_indices(d.d01s, d.d01s.onM));
clear cell_list_map;

% Run with low thresholds then manually check :(
% cell_list_map = map_ei(d.d01s, d.d04, 'master_cell_type', {3}, 'electrode_threshold', 1, 'significant_electrodes', 3, 'corr_threshold', 0.8);
% d.d04.onM_fd01s = cell_list_map(get_cell_indices(d.d01s, d.d01s.onM));
% clear cell_list_map;
d.d04.onM_fd01s = {[] [] 2266 [] [] 3948 4341 [] [] [] 6409 [] [] [] 7144};

cell_list_map = map_ei(d.d01s, d.d05, 'master_cell_type', {3});
d.d05.onM_fd01s = cell_list_map(get_cell_indices(d.d01s, d.d01s.onM));
clear cell_list_map;

d.d04 = read_stim_lisp_output(d.d04);
d.d04.stimulus = parse_stim_rgbs(d.d04.stimulus);

conerun = d.d01s;
conerun.wwrgcs = conerun.onM;

stablerun = d.d05;
stablerun.wwrgcs = stablerun.onM_fd01s;

ud04.wwrun = d.d04;
ud04.wwrun.wwrgcs = ud04.wwrun.onM_fd01s;
ud04.conerun = conerun;
ud04.stablerun = stablerun;
ww_compound_plot(ud04, 'all', 'triggers', ud04.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear conerun stablerun udprintpath


% For J Freeman subunits analysis
save_udraster_data(ud04, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));



%% data004cm U/D onM
udprintpath = [];%printpath('ud', piece);

d.d01_04_05_n = load_data([piece '/d01-04-05-norefit'], loadopts);
d.d01cm = load_data([piece '/d01-04-05-norefit/data001/data001'], loadopts);
d.d04cm = load_data([piece '/d01-04-05-norefit/data004/data004'], loadopts);
d.d05cm = load_data([piece '/d01-04-05-norefit/data005/data005'], staopts);

cell_list_map = map_ei(d.d01s, d.d01_04_05_n, 'master_cell_type', {3});
d.d01_04_05_n.onM_fd01s = cell_list_map(get_cell_indices(d.d01s, d.d01s.onM));
clear cell_list_map;

d.d01cm.onM_fd01s = d.d01_04_05_n.onM_fd01s;
d.d04cm.onM_fd01s = d.d01_04_05_n.onM_fd01s;
d.d05cm.onM_fd01s = d.d01_04_05_n.onM_fd01s;

d.d04cm = read_stim_lisp_output(d.d04cm);
d.d04cm.stimulus = parse_stim_rgbs(d.d04cm.stimulus);

conerun = d.d01s;
conerun.wwrgcs = conerun.onM;

stablerun = d.d05cm;
stablerun.wwrgcs = stablerun.onM_fd01s;

ud04.wwrun = d.d04cm;
ud04.wwrun.wwrgcs = ud04.wwrun.onM_fd01s;
ud04.conerun = conerun;
ud04.stablerun = stablerun;
ww_compound_plot(ud04, 'all', 'triggers', ud04.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear conerun stablerun udprintpath

% For J Freeman subunits analysis
save_udraster_data(ud04, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


%% data004cm C/R onM
crprintpath = [];%printpath('cr', piece);

if ~isfield(d, 'd04cm'), d.d04cm = load_data([piece '/d01-04-05-norefit/data004/data004'], loadopts); end
if ~isfield(d, 'd05cm'), d.d05cm = load_data([piece '/d01-04-05-norefit/data005/data005'], loadopts); end

if ~isfield(d.d04cm, 'onM_fd01s')
    if ~isfield(d, 'd01_04_05'), d.d01_04_05 = load_data([piece '/d01-04-05-norefit'], loadopts); end
    cell_list_map = map_ei(d.d01s, d.d01_04_05, 'master_cell_type', {3});
    d.d04cm.onM_fd01s = cell_list_map(get_cell_indices(d.d01s, d.d01s.onM));
    clear cell_list_map;
end

conerun = d.d01s;
conerun.rgcs = conerun.onM;

rasterrun = d.d04cm;
rasterrun.rgcs = rasterrun.onM_fd01s;

stablerun = d.d05cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun crs crsx ps resnorms residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d04cm d.d01s d.d05cm] = deal(cleared{:});


%% data003 is the same regions as data004 but with opposite (i.e. OFF) polarity; should check those too?

%% data004 searching
if ~isfield(d, 'd04cm'), d.d04cm = load_data([piece '/d01-04-05-norefit/data004/data004'], loadopts); end
map_d04cm_d01s = map_ei(d.d01s, d.d04cm);
d.d04cm = read_stim_lisp_output(d.d04cm);
d.d04cm.stimulus = parse_stim_rgbs(d.d04cm.stimulus);
urgb = d.d04cm.stimulus.urgb;

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d01s;
rasterrun = d.d04cm;
map = map_d04cm_d01s;

% Search all ON Parasols
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
% Interesting ones in here, especially from surround decrement stimulation
% and U/D.  Only one responds in "single cone" stimuli and it has more than
% one set of regions near center, so not really usable.

% OFF Parasols
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
% ON midget 7292
%   OFF parasols 6993, 7666
%
% OFF midget 6796
%   OFF parasols 6798
%
% OFF midget 1741
%   OFF parasols 1591
%
% OFF midgets 2058, 1789, 2117
%   OFF parasols 1906
%
% OFF midget 6527
%   OFF parasol 5506 (not isolated)
%
% OFF parasols 5973, 4577, 2179

% ON Midgets
celltype = 3;
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
% Nothing but the targeted ON midgets

% OFF Midgets
celltype = 4;
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
% ON midget 2795
%   OFF midget 1204 (only 1 cone)
%
% ON midget 7292 & OFF parasols 6993, 7666
%   OFF midget 7516
%
% OFF midgets 1741, 1789, 2058, 2117, 6527, 6796


%% data005 analyze static nonlinearities
d.d05 = load_data([piece '/data005'], staopts);

cell_list_map = map_ei(d.d05s, d.d05, 'master_cell_type', 'all');
d.d05.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d05.offP_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offP));
clear cell_list_map;

d.d05 = load_java_movie(d.d05);
%d.d05 = get_snls_old(d.d05, 1606, 'marks', 'sigstix');
d.d05 = get_snls(d.d05, 1606, 'marks', 'sigstix', 'new', true);
%d.d05 = get_snls(d.d05, 1606, 'marks', 'simple');
plot_snl_(d.d05.stas.snls{74}.gen_signal, d.d05.stas.snls{74}.spikes, 'foa', 0, 'fit', d.d05.stas.snls{74}.fit_params)


% Test RGB
d.d00 = load_data([piece '/data000'], staopts);

cell_list_map = map_ei(d.d05s, d.d00, 'master_cell_type', 'all');
d.d00.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d00.offP_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offP));
clear cell_list_map;

d.d00 = load_java_movie(d.d00);
%d.d00 = get_snls(d.d00, 1727, 'marks', 'sigstix', 'new', true);
d.d00 = get_snls_old(d.d00, 1727, 'marks', 'sigstix', 'new', true);
plot_snl_(d.d00.stas.snls{75}.gen_signal, d.d00.stas.snls{75}.spikes, 'foa', 0, 'fit', d.d00.stas.snls{75}.fit_params)


% For profiling
tic
profile on;
d.d05 = get_snls(d.d05, cell2mat(d.d05.offM_fd05s), 'marks', 'simple', 'new', true);
%d.d05 = get_snls_old(d.d05, cell2mat(d.d05.offM_fd05s), 'marks', 'simple', 'new', true);
profile off;
toc
profile viewer;
plot_snl_(d.d05.stas.snls{74}.gen_signal, d.d05.stas.snls{74}.spikes, 'foa', 0, 'fit', d.d05.stas.snls{74}.fit_params)


%% data007 U/D + C/R offM and offP
udprintpath = [];%printpath('ud', piece);
crprintpath = [];%printpath('cr', piece);

d.d05 = load_data([piece '/data005'], loadopts);
d.d07 = load_data([piece '/data007'], loadopts);
d.d09 = load_data([piece '/data009'], loadopts);

cell_list_map = map_ei(d.d05s, d.d05, 'master_cell_type', 'all');
d.d05.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d05.offP_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offP));
clear cell_list_map;

cell_list_map = map_ei(d.d05s, d.d07, 'master_cell_type', 'all');
d.d07.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d07.offP_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offP));
clear cell_list_map;

cell_list_map = map_ei(d.d05s, d.d09, 'master_cell_type', 'all');
d.d09.offM_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offM));
d.d09.offP_fd05s = cell_list_map(get_cell_indices(d.d05s, d.d05s.offP));
clear cell_list_map;

d.d07 = read_stim_lisp_output(d.d07);
d.d07.stimulus = parse_stim_rgbs(d.d07.stimulus);

conerun = d.d05s;
conerun.wwrgcs = [conerun.offM conerun.offP];

stablerun = d.d09;
stablerun.wwrgcs = [stablerun.offM_fd05s stablerun.offP_fd05s];

ud07.wwrun = d.d07;
ud07.wwrun.wwrgcs = [ud07.wwrun.offM_fd05s ud07.wwrun.offP_fd05s];
ud07.conerun = conerun;
ud07.stablerun = stablerun;
ww_compound_plot(ud07, 'all', 'triggers', ud07.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear conerun stablerun udprintpath


% For J Freeman subunits analysis
save_udraster_data(ud07, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


conerun = d.d05;
conerun.rgcs = [conerun.offM_fd05s conerun.offP_fd05s];
rasterrun = d.d07;
rasterrun.rgcs = [rasterrun.offM_fd05s rasterrun.offP_fd05s];
stablerun = d.d09;
stablerun.rgcs = [stablerun.offM_fd05s stablerun.offP_fd05s];
[rasterrun conerun stablerun crs crsx ps resnorms residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', [6 16], 'printpath', crprintpath, 'scaled_up', 10);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d05s d.d07 d.d09] = deal(cleared{:});

% 1cone paper M- (index 6)
x = 73.5; y = 55.5; del = 29;
axis([x x+del y y+del]);
% 1cone paper P- (index 16)

% Searched finds centers
cell_list_map_d07fd05s = map_ei(d.d05s, d.d07);
cell_list_map_d09fd05s = map_ei(d.d05s, d.d09);
conerun.rgcs = [2506, 2896, 6841, 1981, 2117];
rasterrun.rgcs = cell_list_map_d07fd05s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = cell_list_map_d09fd05s(get_cell_indices(conerun, conerun.rgcs));
[rasterrun conerun stablerun crs crsx ps resnorms residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);

% Searched finds surrounds
conerun.rgcs = [2461, 3331, 1951, 47, 151, 3136, 1786, 1818, 1953, 2056, 2131];
rasterrun.rgcs = cell_list_map_d07fd05s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = cell_list_map_d09fd05s(get_cell_indices(conerun, conerun.rgcs));
[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath, 'surround', true);


%% Unscaled C/Rs for Neuron reviewer 2
subplot(1,4,2)
plot(crsx(1,:)./3.*200, crs(:,:,6)', '.-', 'MarkerSize', 15);
title('OFF Midget')
xlabel('% contrast')

subplot(1,4,4)
plot(crsx(1,:)./3.*200, crs(:,:,16)', '.-', 'MarkerSize', 15);
title('OFF Parasol')
xlabel('% contrast')


%% data007 searching
if ~isfield(d, 'd07'), d.d07 = load_data([piece '/data007'], loadopts); end

cell_list_map_d05s = map_ei(d.d05s, d.d07);
d.d07.all_fd05s = cell_list_map_d05s(get_cell_indices(d.d05s, 'all'));

d.d07 = read_stim_lisp_output(d.d07);
d.d07.stimulus = parse_stim_rgbs(d.d07.stimulus);
urgb = d.d07.stimulus.urgb;

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d05s;
rasterrun = d.d07;

% Nothing for ON Parasols
% Nothing for ON Midgets

% Search all OFF Parasols
celltype = 2;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Center: 2506, 2896, 6841 (d05s)
% Center shared: 1981, 2117 (d05s)
% Surround: 2461, 3331 (d05s)

% Search all OFF Midgets
celltype = 4;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Surround: 1951 (d05s)
% Surround dominated: 47, 151, 3136 (d05s)
% Surround shared: 1786, 1818, 1953, 2056, 2131 (d05s)

% Search all OFF large
celltype = 10;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
end


%% Null stimulus data008
if ~isfield(d, 'd05'), d.d05 = load_data([piece '/data005'], staopts); end
if ~isfield(d, 'd08'), d.d08 = load_data([piece '/data008'], loadopts); end
if ~isfield(d, 'd09'), d.d09 = load_data([piece '/data009'], staopts); end

% Have to be careful not to try to parse polylines because the 2nd stim region
% isn't well formed and this will crash.
d.d08 = read_stim_lisp_output(d.d08, [], false, false);
d.d08.stimulus = parse_stim_rgbs(d.d08.stimulus);

d.d08.mapd05s = map_ei(d.d05s, d.d08);
d.d08.mapd05 = map_ei(d.d05, d.d08);
d.d09.mapd05 = map_ei(d.d05, d.d09);

ploturgbs = {3 2 1 10 11 12; 6 5 4 7 8 9};
titles = {-0.48 -0.24 -0.12 0.12 0.24 0.48};
conerun = d.d05;
rasterrun = d.d08;
map = d.d08.mapd05;
stablerun = d.d09;
stablemap = d.d09.mapd05;

celltype = 1;
% Many of these have nice null rasters, but the RFs look like shit and the
% cone map doesn't match them well...
%
% Maybe is kosher to use only first half of trials, on assumption that
% movement messes things up more in second half; this would improve several
% of them considerably.
%
%  246, 2044, 2343, 7325 best (d05)

celltype = 2;
%  1097 interesting; only surround (d05)
%  376, 1096, 1981, 6031, 7186 best (d05)

celltype = 3;
% Not many
%  1956 best (d05)

celltype = 4;
%  1533, 2026, 2386, 6241 best (d05)

celltype = 5;
%  1321, 2416, 6287 (d05)

celltype = 8;
%  1907 (d05)

for cellid = conerun.cell_types{celltype}.cell_ids(21:40)
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