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
d10.stimulus = convert_stimulus_to_combined_maps(d10.stimulus);
d10.stimulus = parse_stim_rgbs(d10.stimulus);
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
d09cm.stimulus = convert_stimulus_to_combined_maps(d09cm.stimulus);
d09cm.stimulus = parse_stim_rgbs(d09cm.stimulus);
d09cm.stimulus.triggers = d09cm.triggers(1:2:end);
d09cm.mapd06s = map_ei(d06s, d09cm);
d10cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data010-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data010-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], staopts);
d10cm.mapd06s = map_ei(d06s, d10cm);
pieces(piece) = struct('d06s', d06s, 'd09cm', d09cm, 'd10cm', d10cm);
leave(keep_vars{:});

%% DOESN'T WORK!
%% Load 2012-08-21-2
piece = '2012-08-21-2';
d01s = load_data([piece '/streamed/data001/data001'], staopts);
d01s = load_cones_ath(d01s, 2);
d03 = load_data([piece '/data003'], loadopts);
d03 = read_stim_lisp_output(d03);
d03.stimulus = parse_stim_rgbs(d03.stimulus);
d03.stimulus.triggers = d03.triggers(1:2:end);
d03.mapd01s = map_ei(d01s, d03);
d04 = load_data([piece '/data004'], staopts);
d04.mapd01s = map_ei(d01s, d04);
pieces(piece) = struct('d01s', d01s, 'd03', d03, 'd04', d04);
leave(keep_vars{:});

%% Load 2012-09-06-0
piece = '2012-09-06-0';
d04s = load_data([piece, '/streamed/data004/data004'], staopts);
d04s = load_cones_ath(d04s, 'piece');
d08 = load_data([piece '/data008'], loadopts);
d08 = read_stim_lisp_output(d08);
d08.stimulus = parse_stim_rgbs(d08.stimulus);
d08.stimulus.triggers = d08.triggers(1:2:end);
d08.mapd04s = map_ei(d04s, d08);
d09 = load_data([piece '/data009'], staopts);
d09.mapd04s = map_ei(d04s, d09);
pieces(piece) = struct('d04s', d04s, 'd08', d08, 'd09', d09);
leave(keep_vars{:});

%% Load 2012-09-21-2
piece = '2012-09-21-2';
d09 = load_data([piece '/data009'], staopts);
d09 = load_cones_ath(d09, 'piece');
d13 = load_data([piece '/data013'], loadopts);
d13 = read_stim_lisp_output(d13);
d13.stimulus = parse_stim_rgbs(d13.stimulus);
d13.stimulus.triggers = d13.triggers(1:2:end);
d13.mapd09 = map_ei(d09, d13);
d14 = load_data([piece '/data014'], staopts);
d14.mapd09 = map_ei(d09, d14);
pieces(piece) = struct('d09', d09, 'd13', d13, 'd14', d14);
leave(keep_vars{:});





%% FIGURE VERSION; upsampled RF, colors set, etc.
% Allcones, 6 cone off midget 2011-12-13-2 id5162
piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = 5162;
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 6;

localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
localmap = localmask.*rasterrun.stimulus.mapims{1};
[ax rfweights crweights f gof p crs crsx] = allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun, 'rfopts', {'scaled_up', 10, 'colors', zeros(6,3)}, 'xscaleopts', {'colors', zeros(6,3), 'Markers', 'dxs^o+', 'MarkerSize', 5}, 'marksthresh', 3.5);

% Unscaled C/R for Neuron reviewer 2
figure()
plot(crsx(1,:)./3.*200, crs(end-5:end,:)', '.-', 'MarkerSize', 15);
title('OFF Midget')
xlabel('% contrast')
ylabel('# spikes')

leave(keep_vars{:});




%% Allcones, 6 cone off midget 2011-12-13-2 id5162
piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = 5162;
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 6;

localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
localmap = localmask.*rasterrun.stimulus.mapims{1};
allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);

leave(keep_vars{:});


%% Allcone, 16 (11) cone off midget 2011-12-13-2 id1351
piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = 1351;
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 6;

localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
localmap = localmask.*rasterrun.stimulus.mapims{1};
localmap(localmap == 10) = 0;
allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);

leave(keep_vars{:});


%% Allcone, 8 (7) cone off midget 2011-12-13-2 id2586
% Needs to have one region removed
piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = 3586;
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 6;

localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
localmap = localmask.*rasterrun.stimulus.mapims{1};
localmap(localmap == 13) = 0;
allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);

leave(keep_vars{:});


%% Allcones, 6 cone off midget 2012-08-21-2 id887
% No stability as yet, doesn't look very good in new analysis 5162 nicer

% piece = pieces('2012-08-21-2');
% conerun = piece.d01s;
% conerun.rgcs = 887;
% rasterrun = piece.d03;
% stablerun = piece.d04;
% rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
% stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
% padfactor = 5;
% 
% localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
% localmap = localmask.*rasterrun.stimulus.mapims{1};
% localmap(localmap == 13) = 0;
% allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);
% 
% leave(keep_vars{:});


%% Allcones, 13 (10) cone off midget 2012-09-21-2 id5086
% Stability not great
piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = 5086;
rasterrun = piece.d13;
stablerun = piece.d14;
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
padfactor = 5;

localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
localmap = localmask.*rasterrun.stimulus.mapims{1};
allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);
leave(keep_vars{:});


%% Allcones, 18 (14) cone off midget 2012-09-21-2 id6346
% Glaring missing cone, doesn't look too impressive in new analysis

% piece = pieces('2012-09-21-2');
% conerun = piece.d09;
% conerun.rgcs = 6346;
% rasterrun = piece.d13;
% stablerun = piece.d14;
% rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
% stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
% padfactor = 4.2;
% ar = 1.25;
% 
% localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor, 'az_aspect_ratio', ar);
% localmap = localmask.*rasterrun.stimulus.mapims{1};
% localmap(localmap == 16) = 0;
% allcones_plot(conerun, rasterrun, localmap, padfactor, 'ar', ar, 'stablerun', stablerun);
% 
% leave(keep_vars{:});


%% Allcones, 18 (14) cone off midget 2012-09-21-2 id6406
% Pretty good
piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = 6406;
rasterrun = piece.d13;
stablerun = piece.d14;
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
padfactor = 7.5;

localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
localmap = localmask.*rasterrun.stimulus.mapims{1};
localmap(localmap == 15) = 0;
allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);
leave(keep_vars{:});


%% Allcones, 4 cone off midget 2012-08-21-2 id6005
% 7682 nicer

% piece = pieces('2012-08-21-2');
% conerun = piece.d01s;
% conerun.rgcs = 6005;
% rasterrun = piece.d03;
% stablerun = piece.d04;
% rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
% stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
% padfactor = 6;
% 
% localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
% localmap = localmask.*rasterrun.stimulus.mapims{1};
% allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);

% leave(keep_vars{:});


%% Allcones, 4 cone off midget 2012-08-21-2 id7682
piece = pieces('2012-08-21-2');
conerun = piece.d01s;
conerun.rgcs = 7682;
rasterrun = piece.d03;
stablerun = piece.d04;
rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 1.8;

localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
localmap = localmask.*rasterrun.stimulus.mapims{1};
allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);

leave(keep_vars{:});


%% Allcones, 5 cone off midget 2012-09-06-0 id5746


%% Allcones, 13 (5) cone off midget 2012-09-06-0 id436
% Stability issues

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

%% Picking regions, checking fits, etc.
piece = pieces('2012-09-06-0');
conerun = piece.d04s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(14);
rasterrun = piece.d08;
stablerun = piece.d09;
rasterrun.rgcs = rasterrun.mapd04s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd04s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 4;
ar = 1.4;
badregions = [];

% Get locally masked off stimulus map
localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor, 'az_aspect_ratio', ar);
localmap = localmask.*rasterrun.stimulus.mapims{1};

% Further crop bad regions?
for crop = badregions
    localmap(localmap == crop) = 0;
end

% Get cones stimulated
cbi = recover_cones_stimulated(conerun, localmap)

% Check for proper regions
localpoly = nyc2poly(map2manhattan(localmap));
ax = plot_rf_stimmap(conerun, conerun.rgcs, localpoly, 'fit', false, 'az_pad_factor', padfactor, 'az_aspect_ratio', ar);
bounds = axis(ax);
plot_cone_mosaic(conerun, 'cone_roi', vertcat(cbi{:}), 'label', true, 'fig_or_axes', ax, 'clear', false, 'roi_highlight', false, 'bounds', bounds);
axis(bounds);

%%
[ax rfweights crweights f gof] = allcones_plot(conerun, rasterrun, localmap, padfactor, 'ar', ar, 'stablerun', stablerun);

leave(keep_vars{:});


%% Scatterplot - gather data
sepfigs = false;
if ~sepfigs, figure; end

startind = 0;
piece = pieces('2011-12-13-2');
rgcs        = [1321 1351  2251 3586          4576 5162];
padfactors  = [   6    6     6    6             6    6];
badregionss = {  [] [10] [8 9] [13] [6 7 8 10 11]   []};
conerun = piece.d08s;
rasterrun = piece.d10;
stablerun = piece.d14s;
for i = 1:length(rgcs)
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    badregions = badregionss{i};
    rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
    stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
    
    localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
    localmap = localmask.*rasterrun.stimulus.mapims{1};
    for r = badregions, localmap(localmap == r) = 0; end
    if sepfigs, figure; end
    [axs rfweights{i+startind} crweights{i+startind} f{i+startind} gof(i+startind)] = allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);        
    drawnow
end

% 2536 xscale fitting not great
startind = 6;
piece = pieces('2012-04-13-1');
rgcs       = [2536 6136];
padfactors = [   5    5];
conerun = piece.d06s;
conerun.rgcs = 2536;
rasterrun = piece.d09cm;
stablerun = piece.d10cm;
for i = length(rgcs)
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    stablerun.rgcs = stablerun.mapd06s(get_cell_indices(conerun, conerun.rgcs));
    rasterrun.rgcs = stablerun.rgcs;
    
    localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
    localmap = localmask.*rasterrun.stimulus.mapims{1};
    if sepfigs, figure; end
    [axs rfweights{i+startind} crweights{i+startind} f{i+startind} gof(i+startind)] = allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);        
    drawnow
end

% 887 could be cut for stability...
startind = 8;
piece = pieces('2012-08-21-2');
rgcs       = [391 887 2266 2828 5447 7682];
padfactors = [5.8   5    5    4    5  1.8];
conerun = piece.d01s;
conerun.rgcs = 7682;
rasterrun = piece.d03;
stablerun = piece.d04;
for i = [1 3:length(rgcs)]
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
    stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));

    localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
    localmap = localmask.*rasterrun.stimulus.mapims{1};
    if sepfigs, figure; end
    [axs rfweights{i+startind} crweights{i+startind} f{i+startind} gof(i+startind)] = allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);        
    drawnow
end

% 904 stability issues
% 7022 minor stability issues but nice
startind = 14;
piece = pieces('2012-09-06-0');
rgcs       = [904 5103 5746 6031 6317 7022];
padfactors = [  6    4    4    5    4    7];
conerun = piece.d04s;
conerun.rgcs = 7022;
rasterrun = piece.d08;
stablerun = piece.d09;
for i = 2:length(rgcs)
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    rasterrun.rgcs = rasterrun.mapd04s(get_cell_indices(conerun, conerun.rgcs));
    stablerun.rgcs = stablerun.mapd04s(get_cell_indices(conerun, conerun.rgcs));

    localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
    localmap = localmask.*rasterrun.stimulus.mapims{1};
    if sepfigs, figure; end
    [axs rfweights{i+startind} crweights{i+startind} f{i+startind} gof(i+startind)] = allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun);        
    drawnow
end

% 5086 some stability issues but nice
% 5491 stability a bit suspect
% 7006 xscale fit not great
startind = 20;
piece = pieces('2012-09-21-2');
rgcs        = [5086 5491 6346 6406 7006 7696];
padfactors  = [   5    7  4.2    8    4    5];
badregionss = {  []   [] [16] [15]   []   []};
ars         = [   1    1  1.2    1    1    1];
conerun = piece.d09;
rasterrun = piece.d13;
stablerun = piece.d14;
urgbs = {[1 1 1].*-0.48};
for i = [1 3:4 6]
    conerun.rgcs = rgcs(i);
    padfactor = padfactors(i);
    ar = ars(i);
    badregions = badregionss{i};
    rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
    stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));

    localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor, 'az_aspect_ratio', ar);
    localmap = localmask.*rasterrun.stimulus.mapims{1};
    for r = badregions, localmap(localmap == r) = 0; end
    if sepfigs, figure; end
    [axs rfweights{i+startind} crweights{i+startind} f{i+startind} gof(i+startind)] = allcones_plot(conerun, rasterrun, localmap, padfactor, 'stablerun', stablerun, 'ar', ar);
    drawnow
end


%%
save('/home/peter/Dropbox/1cone paper/draftfigs3/allcones/scatterweights21b', 'crweights', 'rfweights', 'f', 'gof');


%% Scatterplot
load('/home/peter/Dropbox/1cone paper/draftfigs3/allcones/scatterweights21b');

figure;
for i = 1:length(crweights)
    if isempty(f{i}), continue; end
    normcrweights{i} = f{i}.a * crweights{i};
end

cla(); hold on;

% Cells shown
h = plot(rfweights{2},  normcrweights{2},  'ko');
h = plot(rfweights{4},  normcrweights{4},  'ks');
h = plot(rfweights{6},  normcrweights{6},  'kd');
h = plot(rfweights{21}, normcrweights{21}, 'k^');
h = plot(rfweights{24}, normcrweights{24}, 'kv');

% The rest
for i = [1 3 5 7:20 22 23 25:length(rfweights)]
    if isempty(rfweights{i}), continue; end
    h = plot(rfweights{i}, normcrweights{i}, 'k.');
end

plot([0 1], [0 1], '--k');

axis equal
axis([-0.1 1.1 0 1.1]);