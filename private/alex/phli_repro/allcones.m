%% Setup
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');
keep_vars = {'pieces'; 'loadopts'; 'staopts'; 'keep_vars'};

%% Load 2011-12-13-2
piece = '2011-12-13-2';
d08s = load_data([piece '/streamed/data008-0/data008-0'], staopts);
d08s = load_cones_ath(d08s, 'bayes');
d10 = load_data([piece '/data010'], loadopts);
d10 = read_stim_lisp_output(d10, ':2011-12-13-2_f08_allcones');
d10.stimulus.triggers = d10.triggers(1:2:end);
d10.mapd08s = map_ei(d08s, d10);
d14s = load_data([piece '/streamed/data014-0/data014-0'], staopts);
d14s.mapd08s = map_ei(d08s, d14s);
pieces(piece) = struct('d08s', d08s, 'd10', d10, 'd14s', d14s);
leave(keep_vars{:});

%% Load 2012-04-13-1
piece = '2012-04-13-1';
d06s = load_data([piece '/streamed/data006/data006'], staopts);
d06s = load_cones_ath(d06s, piece);
d09cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data009-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data009-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
d09cm = read_stim_lisp_output(d09cm, ':2012-04-13-1:f06_allcones');
d09cm.stimulus.triggers = d09cm.triggers(1:2:end);
d09cm.mapd06s = map_ei(d06s, d09cm);
d10cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data010-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data010-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], staopts);
d10cm.mapd06s = map_ei(d06s, d10cm);
pieces(piece) = struct('d06s', d06s, 'd09cm', d09cm, 'd10cm', d10cm);
leave(keep_vars{:});

%% Load 2012-08-21-2
piece = '2012-08-21-2';
d01s = load_data([piece '/streamed/data001/data001'], staopts);
d01s = load_cones_ath(d01s, 2);
d03 = load_data([piece '/data003'], loadopts);
d03 = read_stim_lisp_output(d03);
d03.stimulus.triggers = d03.triggers(1:2:end);
d03.mapd01s = map_ei(d01s, d03);
d04 = load_data([piece '/data004'], staopts);
d04.mapd01s = map_ei(d01s, d04);
pieces(piece) = struct('d01s', d01s, 'd03', d03, 'd04', d04);
leave(keep_vars{:});

%% Load 2012-09-06-0
piece = '2012-09-06-0';
d04s = load_data([piece, '/streamed/data004/data004'], staopts);
d04s = load_cones_ath(d04s,1);
d08 = load_data([piece '/data008'], loadopts);
d08 = read_stim_lisp_output(d08);
d08.stimulus.triggers = d08.triggers(1:2:end);
d08.mapd04s = map_ei(d04s, d08);
d09 = load_data([piece '/data009'], staopts);
d09.mapd04s = map_ei(d04s, d09);
pieces(piece) = struct('d04s', d04s, 'd08', d08, 'd09', d09);
leave(keep_vars{:});

%% Load 2012-09-21-2
piece = '2012-09-21-2';
d09 = load_data([piece '/data009'], staopts);
d09 = load_cones_ath(d09,1);
d13 = load_data([piece '/data013'], loadopts);
d13 = read_stim_lisp_output(d13);
d13.stimulus.triggers = d13.triggers(1:2:end);
d13.mapd09 = map_ei(d09, d13);
d14 = load_data([piece '/data014'], staopts);
d14.mapd09 = map_ei(d09, d14);
pieces(piece) = struct('d09', d09, 'd13', d13, 'd14', d14);
leave(keep_vars{:});


%% Allcones, 6 cone off midget 2011-12-13-2 id5162
piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = 5162;
rasterrun = piece.d10;
stablerun = piece.d14s;
urgbs = {[1 1 1].*-0.48};

rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
allcones_plot(rasterrun, conerun, 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
leave(keep_vars{:});


%% Allcone, 16 (11) cone off midget 2011-12-13-2 id1351
% piece = pieces('2011-12-13-2');
% conerun = piece.d08s;
% conerun.rgcs = 1351;
% rasterrun = piece.d10;
% stablerun = piece.d14s;
% urgbs = {[1 1 1].*-0.48};
% 
% rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
% stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
% allcones_plot(rasterrun, conerun, 'az_pad_factor', 6, 'mapweightindices', [1:8 10:17], 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
% leave(keep_vars{:});


%% Allcone, 8 (7) cone off midget 2011-12-13-2 id2586
% Needs to have one region removed
piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = 3586;
rasterrun = piece.d10;
stablerun = piece.d14s;
urgbs = {[1 1 1].*-0.48};

rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
allcones_plot(rasterrun, conerun, 'mapweightindices', [1:7 9], 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
leave(keep_vars{:});


%% Allcones, 6 cone off midget 2012-08-21-2 id887
% No stability as yet, 5162 nicer
% piece = pieces('2012-08-21-2');
% conerun = piece.d01s;
% conerun.rgcs = 887;
% rasterrun = piece.d03;
% stablerun = piece.d04;
% urgbs = {[1 1 1].*-0.48};
% 
% rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
% stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
% allcones_plot(rasterrun, conerun, 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10}, 'az_pad_factor', 5);
% leave(keep_vars{:});


%% Allcones, 13 (10) cone off midget 2012-09-21-2 id5086
piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = 5086;
rasterrun = piece.d13;
stablerun = piece.d14;
urgbs = {[1 1 1].*-0.48};

rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
allcones_plot(rasterrun, conerun, 'az_pad_factor', 5, 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
leave(keep_vars{:});


%% Allcones, 18 (14) cone off midget 2012-09-21-2 id6346
% Nice, but glaring missing cone
% piece = pieces('2012-09-21-2');
% conerun = piece.d09;
% conerun.rgcs = 6346;
% rasterrun = piece.d13;
% stablerun = piece.d14;
% urgbs = {[1 1 1].*-0.48};
% 
% rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
% stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
% allcones_plot(rasterrun, conerun, 'az_pad_factor', 4.2, 'az_aspect_ratio', 1.25, 'mapweightindices', [1:15 17:19], 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
% leave(keep_vars{:});


%% Allcones, 18 (14) cone off midget 2012-09-21-2 id6406
% Pretty good
piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = 6406;
rasterrun = piece.d13;
stablerun = piece.d14;
urgbs = {[1 1 1].*-0.48};

rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
allcones_plot(rasterrun, conerun, 'az_pad_factor', 7.5, 'mapweightindices', [1:12 14:16 18:20], 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
leave(keep_vars{:});


%% Allcones, 4 cone off midget 2012-08-21-2 id6005
% 7682 nicer
% piece = pieces('2012-08-21-2');
% conerun = piece.d01s;
% conerun.rgcs = 6005;
% rasterrun = piece.d03;
% stablerun = piece.d04;
% urgbs = {[1 1 1].*-0.48};
% 
% rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
% stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
% allcones_plot(rasterrun, conerun, 'az_pad_factor', 6, 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
% leave(keep_vars{:});


%% Allcones, 4 cone off midget 2012-08-21-2 id7682
piece = pieces('2012-08-21-2');
conerun = piece.d01s;
conerun.rgcs = 7682;
rasterrun = piece.d03;
stablerun = piece.d04;
urgbs = {[1 1 1].*-0.48};

rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
allcones_plot(rasterrun, conerun, 'az_pad_factor', 1.8, 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
leave(keep_vars{:});


%% Allcones, 5 cone off midget 2012-09-06-0 id5746


%% Allcones, 13 (5) cone off midget 2012-09-06-0 id436
% piece = pieces('2012-09-06-0');
% conerun = piece.d04s;
% conerun.rgcs = 436;
% rasterrun = piece.d08;
% stablerun = piece.d09;
% urgbs = {[1 1 1].*-0.48};
% 
% rasterrun.rgcs = rasterrun.mapd04s(get_cell_indices(conerun, conerun.rgcs));
% stablerun.rgcs = stablerun.mapd04s(get_cell_indices(conerun, conerun.rgcs));
% allcones_plot(rasterrun, conerun, 'urgbs', urgbs, 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
% leave(keep_vars{:});


%% Resize figure widths to get same raster axis width
desired_aspect_ratio = 0.75;

fig = gcf();
c = get(fig, 'Children');
rastax = c(2); % Just grab the most convenient one

% Get currrent position in absolute units
u = get(rastax, 'Units');
set(rastax, 'Units', 'inch');
p = get(rastax, 'Position');
set(rastax, 'Units', u);

% Calculate adjustment
ar = p(3)./p(4);
xscale = desired_aspect_ratio ./ ar;
p = get(fig, 'Position');
p(3) = p(3) .* xscale;
p(3:4) = p(3:4) ./ 3; % Shrink so stupid shit will fit in printout!!
set(fig, 'Position', p);


%% Scatterplot - gather data
urgbs = {[1 1 1].*-0.48};
opts = {'urgbs', urgbs, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10}};
load('/Users/peterli/Dropbox/1cone paper/draftfigs/fig4_allcones/scatterweights');

figs = [];
startind = 0;
piece = pieces('2011-12-13-2');
rgcs       = [1321        1351         2251    3586              4576 5162];
padfactors = [   6           6            6       6                 6    6];
mapweights = {  [] [1:8 10:17] [1:12 14:19] [1:7 9] [1:8 10 12 13 15]   []};
conerun = piece.d08s;
rasterrun = piece.d10;
stablerun = piece.d14s;
for i = 1:length(rgcs)
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    mapweightindices = mapweights{i};
    rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
    stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
    moreopts = {'az_pad_factor', padfactor, 'stablerun', stablerun};
    if ~isempty(mapweightindices), moreopts = [moreopts 'mapweightindices' mapweightindices]; end
    [figs(i) rfweights(i+startind) rasterweights(i+startind)] = allcones_plot(rasterrun, conerun, opts{:}, moreopts{:});
end
close(figs);

figs = [];
startind = 6;
piece = pieces('2012-04-13-1');
rgcs = [2536 6136];
padfactors = [5 5];
conerun = piece.d06s;
conerun.rgcs = 2536;
rasterrun = piece.d09cm;
stablerun = piece.d10cm;
for i = 1:length(rgcs)
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    stablerun.rgcs = stablerun.mapd06s(get_cell_indices(conerun, conerun.rgcs));
    rasterrun.rgcs = stablerun.rgcs;
    moreopts = {'az_pad_factor', padfactor, 'stablerun', stablerun};
    [figs(i) rfweights(i+startind) rasterweights(i+startind)] = allcones_plot(rasterrun, conerun, opts{:}, moreopts{:});
end
close(figs);

figs = [];
startind = 8;
piece = pieces('2012-08-21-2');
rgcs       = [391 887 2266 2828 5447 7682];
padfactors = [5.8   5    5    4    5  1.8];
conerun = piece.d01s;
conerun.rgcs = 7682;
rasterrun = piece.d03;
stablerun = piece.d04;
for i = 1:length(rgcs)
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
    stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
    moreopts = {'az_pad_factor', padfactor, 'stablerun', stablerun};
    [figs(i) rfweights(i+startind) rasterweights(i+startind)] = allcones_plot(rasterrun, conerun, opts{:}, moreopts{:});
end
close(figs);

figs = [];
startind = 14;
piece = pieces('2012-09-06-0');
rgcs       = [904 5103 5746 6031 6317 7022];
padfactors = [  6    4    4    5    4    7];
conerun = piece.d04s;
conerun.rgcs = 7022;
rasterrun = piece.d08;
stablerun = piece.d09;
for i = 1:length(rgcs)
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    rasterrun.rgcs = rasterrun.mapd04s(get_cell_indices(conerun, conerun.rgcs));
    stablerun.rgcs = stablerun.mapd04s(get_cell_indices(conerun, conerun.rgcs));
    moreopts = {'az_pad_factor', padfactor, 'stablerun', stablerun};
    [figs(i) rfweights(i+startind) rasterweights(i+startind)] = allcones_plot(rasterrun, conerun, opts{:}, moreopts{:});
end
close(figs);

figs = [];
startind = 20;
piece = pieces('2012-09-21-2');
rgcs       = [5086 5491         6346               6406 7006 7696];
padfactors = [   5    7          4.2                  8    4    5];
mapweights = {  []   [] [1:15 17:19] [1:12 14:16 18:20]   []   []};
conerun = piece.d09;
rasterrun = piece.d13;
stablerun = piece.d14;
urgbs = {[1 1 1].*-0.48};
for i = [1 2 4:6]
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    mapweightindices = mapweights{i};
    rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
    stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
    moreopts = {'az_pad_factor', padfactor, 'stablerun', stablerun};
    if ~isempty(mapweightindices), moreopts = [moreopts 'mapweightindices' mapweightindices]; end
    [figs(i) rfweights(i+startind) rasterweights(i+startind)] = allcones_plot(rasterrun, conerun, opts{:}, moreopts{:});
end
i = 3;
conerun.rgcs = rgcs(i);
padfactor = padfactors(i);
mapweightindices = mapweights{i};
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
moreopts = {'az_pad_factor', padfactor, 'stablerun', stablerun};
moreopts = [moreopts 'mapweightindices' mapweightindices 'az_aspect_ratio' 1.2];
[figs(i) rfweights(i+startind) rasterweights(i+startind)] = allcones_plot(rasterrun, conerun, opts{:}, moreopts{:});
close(figs);

save('/Users/peterli/Dropbox/1cone paper/draftfigs/fig4_allcones/scatterweights', 'rasterweights', 'rfweights');


%% Scatterplot - very basic analysis, independently normalized
% load('/Users/peterli/Dropbox/1cone paper/draftfigs/fig4_allcones/scatterweights');
% 
% normrfweights     = rfweights;
% normrasterweights = rasterweights;
% for i = 1:length(normrfweights)
%     normrfweights{i} = normrfweights{i} - min(rfweights{i});
%     normrfweights{i} = normrfweights{i} ./ max(rfweights{i});
%     normrasterweights{i} = normrasterweights{i}  - min(normrasterweights{i});
%     normrasterweights{i} = normrasterweights{i} ./ max(normrasterweights{i});
% end
% 
% colors = jet(length(normrfweights));
% cla(); hold on;
% for i = [1:19 20:26]%:length(normrfweights)
%     h = plot(normrfweights{i}, normrasterweights{i}, 'o');
%     set(h, 'Color', colors(i,:));
% end


%% Scatterplot - very basic analysis, normalized to RF max, BW
% I don't have a good baseline firing rate yet, so for now normalize to
% minimum rfweight at bottom end.
load('/Users/peterli/Dropbox/1cone paper/draftfigs/fig4_allcones/scatterweights');

normrfweights     = rfweights;
normrasterweights = rasterweights;
for i = 1:length(normrfweights)
    [minrfweight, minweighti] = min(normrfweights{i});
    normrfweights{i} = normrfweights{i} - minrfweight;
    [maxrfweight, maxweighti] = max(normrfweights{i});
    normrfweights{i} = normrfweights{i} ./ maxrfweight;

    normrasterweights{i} = normrasterweights{i}  - normrasterweights{i}(minweighti);
    normrasterweights{i} = normrasterweights{i} ./ normrasterweights{i}(maxweighti);
end

cla(); hold on;

% Cells shown
plot(normrfweights{6},  normrasterweights{6},  'ks');
plot(normrfweights{21}, normrasterweights{21}, 'ko');
plot(normrfweights{24}, normrasterweights{24}, 'k^');

% The rest
for i = [1:5 7:20 22 23 25:26]%:length(normrfweights)
    h = plot(normrfweights{i}, normrasterweights{i}, 'k.');
end

plot([0 1], [0 1], '--k');