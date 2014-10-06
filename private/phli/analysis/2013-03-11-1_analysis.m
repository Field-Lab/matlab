% Basics
piece = '2013-03-11-1';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};

% This one has a full OFF amacrine cell mosaic as well as some interesting
% weird cells

%% White noise runs and cell picks

% Movement, but some good ones
d.d06s = load_data(fullfile(piece, 'streamed/data006/data006'), staopts);
d.d06s.offM = [496 511 1126 1234 1666 2191 2341 2956 3107 3334 3391 3796 4141 5236 6871 6886];

d.d06 = load_data(fullfile(piece, 'data006'), staopts);
cell_list_map = map_ei(d.d06s, d.d06, 'master_cell_type', {4});
d.d06.offM = cell_list_map(get_cell_indices(d.d06s, d.d06s.offM));
d.d06.offM = [d.d06.offM{:}];
clear cell_list_map;

jfcones_d06 = [1382, 1397, 1365, 1372
1264, 1166, 1227, 1235
1018, 963, 964, 1046
1358, 1381, 1328, 1389
1045, 1090, 1030, 1149
966, 1004, 972, 930
907, 926, 900, 920
810, 847, 839, 828
603, 631, 578, 620
515, 554, 525, 545
470, 500, 451, 466
220, 252, 251, 212
302, 295, 348, 338
393, 423, 417, 438
690, 691, 657, 623
794, 833, 817, 850];


% Fairly stable
d.d10s = load_data(fullfile(piece, 'streamed/data010/data010'), staopts);
d.d10s.offM = [331 424 511 1069 1232 1666 2191 2356 2881 3106 3333 3391 3991 4141 4876 5343 5626 6076 7066 7547];

d.d10 = load_data(fullfile(piece, 'data010'), staopts);
cell_list_map = map_ei(d.d10s, d.d10, 'master_cell_type', {4});
d.d10.offM = cell_list_map(get_cell_indices(d.d10s, d.d10s.offM));
d.d10.offM = [d.d10.offM{:}];
clear cell_list_map;


% Fairly stable
d.d13s = load_data(fullfile(piece, 'streamed/data013/data013'), staopts);
d.d13s.offM = [331 428 511 1067 1085 1246 1666 2356 2881 3106 3334 3391 3796 4141 4801 5343 6574 7067 7547];

d.d13 = load_data(fullfile(piece, 'data013'), staopts);
cell_list_map = map_ei(d.d13s, d.d13, 'master_cell_type', {4});
d.d13.offM = cell_list_map(get_cell_indices(d.d13s, d.d13s.offM));
d.d13.offM = [d.d13.offM{:}];
clear cell_list_map;

jfcones_d13 = [1003, 1029, 1024, 996
1300, 1313, 1322, 1338
1164, 1156, 1137, 1178
984, 977, 1000, 1026
1312, 1328, 1319, 1327
1101, 1135, 1047, 1035
904, 928, 935, 951
845, 877, 851, 876
618, 654, 635, 665
553, 596, 626, 646
378, 373, 318, 390
282, 304, 329, 275
155, 173, 201, 154
835, 855, 883, 856
1086, 1132, 1098, 1111
1172, 1198, 1199, 1188
272, 273, 231, 238
538, 569, 520, 639
750, 762, 848, 854];


%% data009 U/D
udprintpath = [];%printpath('ud', piece);

d.d09 = load_data([piece '/data009'], loadopts);
cell_list_map = map_ei(d.d06s, d.d09, 'master_cell_type', {4});
d.d09.offM_fd06s = cell_list_map(get_cell_indices(d.d06s, d.d06s.offM));
clear cell_list_map;

cell_list_map = map_ei(d.d06s, d.d10, 'master_cell_type', {4});
d.d10.offM_fd06s = cell_list_map(get_cell_indices(d.d06s, d.d06s.offM));
clear cell_list_map;

d.d09 = read_stim_lisp_output(d.d09);
d.d09.stimulus = parse_stim_rgbs(d.d09.stimulus);

conerun = d.d06;
conerun.wwrgcs = conerun.offM;

stablerun = d.d10;
stablerun.wwrgcs = stablerun.offM_fd06s;

ud09.wwrun = d.d09;
ud09.wwrun.wwrgcs = ud09.wwrun.offM_fd06s;
ud09.conerun = conerun;
ud09.stablerun = stablerun;
ww_compound_plot(ud09, 'all', 'triggers', ud09.wwrun.triggers(1:2:end), 'printpath', udprintpath, 'closefigs', true);
clear conerun stablerun udprintpath


%% data012 searching
if ~isfield(d, 'd12'), d.d12 = load_data([piece '/data012'], loadopts); end
d.d12.mapd10 = map_ei(d.d10, d.d12);
d.d12 = read_stim_lisp_output(d.d12);
d.d12.stimulus = parse_stim_rgbs(d.d12.stimulus);
urgb = d.d12.stimulus.urgb;

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};

conerun = d.d10;
rasterrun = d.d12;
map = d.d12.mapd10;

% Search all OFF amacrine
celltype = 6;
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Nothing :(


%% Null stimulus data011
if ~isfield(d, 'd11'), d.d11 = load_data([piece '/data011'], loadopts); end

% Have to be careful not to try to parse polylines because the 2nd stim region
% isn't well formed and this will crash.
d.d11 = read_stim_lisp_output(d.d11, [], false, false);
d.d11.stimulus = parse_stim_rgbs(d.d11.stimulus);

d.d11.mapd10 = map_ei(d.d10, d.d11);
d.d13.mapd10 = map_ei(d.d10, d.d13);

ploturgbs = {9 8 7 10 11 12; 3 2 1 4 5 6};
titles = {-0.48 -0.24 -0.12 0.12 0.24 0.48};
conerun = d.d10;
rasterrun = d.d11;
map = d.d11.mapd10;
stablerun = d.d13;
stablemap = d.d13.mapd10;

celltype = 1;
% Pretty good

celltype = 2;

celltype = 3;

celltype = 4;
% Missed some good ones, but these are pretty good:
% 1070, 2356, 6076, 7547

celltype = 13;


for cellid = conerun.cell_types{celltype}.cell_ids
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