%% Setup

loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};

pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');

keep_vars = {'pieces'; 'loadopts'; 'staopts'; 'keep_vars'};


%% Load 2012-09-24-5
piece = '2012-09-24-5';
% d01s = load_data([piece '/data001'], streamedopts);
d03s = load_data([piece '/data003'], staopts);
d05s = load_data([piece '/data005'], staopts);
d03cm = load_data([piece '/d03-06-07-norefit/data003/data003'], staopts);
d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], staopts);

d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d03_06_07.mapd03s = map_ei(d03s, d03_06_07);
d03_06_07.mapd05s = map_ei(d05s, d03_06_07);

d06cm = read_stim_lisp_output(d06cm);
d06cm.stimulus = parse_stim_rgbs(d06cm.stimulus);
d06cm.stimulus = parse_cr_rgbs(d06cm.stimulus);

pieces(piece) = struct('d03s', d03s, 'd05s', d05s, 'd03cm', d03cm, 'd06cm', d06cm, 'd07cm', d07cm, 'd03_06_07', d03_06_07);
leave(keep_vars{:});


%% Divergence 1, 2012-09-24-5
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d06cm;
maprun = piece.d03_06_07;
offM = 2851;
onM = 2855;
offP = 2521;
onP = 2988;
cones2plot = 3:4;

bounds = [110 165 70 120];

colors = [1 0 0; 0 0 1];
rasteropts = {'hist_line_width' 2};
rfopts = {'fit', false, 'autozoom', false, 'scaled_up', 10};

urgb = rasterrun.stimulus.urgb;
singleurgbs = find(urgb.singles);
for i = 1:length(cones2plot);
    intensities = urgb.intensities(cones2plot(i),singleurgbs);
    offurgbs(i) = singleurgbs(intensities == -0.48*3);
    onurgbs(i)  = singleurgbs(intensities ==  0.48*3);
end

triggers = rasterrun.triggers(1:2:end);

cellid = maprun.mapd03s{get_cell_indices(conerun, offM)};
urgbs = cell(3,8);
urgbs(3,[2 1]) = {offurgbs(1) offurgbs(2)};
hist_colors = cell(3,8);
hist_colors{3,2}{1} = colors(1,:);
hist_colors{3,1}{1} = colors(2,:);
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});

cellid = maprun.mapd03s{get_cell_indices(conerun, offP)};
urgbs = cell(3,8);
urgbs(3,[4 3]) = {offurgbs(1) offurgbs(2)};
hist_colors = cell(3,8);
hist_colors{3,4}{1} = colors(1,:);
hist_colors{3,3}{1} = colors(2,:);
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});

cellid = maprun.mapd03s{get_cell_indices(conerun, onM)};
urgbs = cell(3,8);
urgbs(3,[6 5]) = {onurgbs(1) onurgbs(2)};
hist_colors = cell(3,8);
hist_colors{3,6}{1} = colors(1,:);
hist_colors{3,5}{1} = colors(2,:);
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});


doubleurgbs = find(urgb.doubles & ~urgb.uds);
intensities = sum(urgb.intensities(cones2plot,doubleurgbs));

cellid = maprun.mapd03s{get_cell_indices(conerun, onP)};
urgbs = cell(3,12);
urgbs(3,[11 10]) = {onurgbs(1) onurgbs(2)};
urgbs{3,12} = doubleurgbs(intensities == 0.48*6);
hist_colors = cell(3,12);
hist_colors{3,11}{1} = colors(1,:);
hist_colors{3,10}{1} = colors(2,:);
hist_colors{3,12}{1} = [1 0 1];
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});

sanesubplot(3,4,{1:2 1});
plot_rf_stimmap(conerun, offM, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);
sanesubplot(3,4,{1:2 2});
plot_rf_stimmap(conerun, offP, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);
sanesubplot(3,4,{1:2 3});
plot_rf_stimmap(conerun, onM, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);
sanesubplot(3,4,{1:2 4});
plot_rf_stimmap(conerun, onP, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);

set(gcf, 'PaperType', 'tabloid');
orient(gcf, 'landscape');

leave(keep_vars{:});


%% Divergence 2, 2012-09-24-5
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d06cm;
maprun = piece.d03_06_07;
offM1 = 3646;
offM2 = 3631;
onM_d05s = 3738;
offP = 3616;
onP = 3739;
cones2plot = 2:3;

bounds = [160 220 15 85];

colors = [1 0 0; 0 0 1];
rasteropts = {'hist_line_width' 2};
rfopts = {'fit', false, 'autozoom', false, 'scaled_up', 10};

urgb = rasterrun.stimulus.urgb;
singleurgbs = find(urgb.singles);
for i = 1:length(cones2plot);
    intensities = urgb.intensities(cones2plot(i),singleurgbs);
    offurgbs(i) = singleurgbs(intensities == -0.48*3);
    onurgbs(i)  = singleurgbs(intensities ==  0.48*3);
end

triggers = rasterrun.triggers(1:2:end);

f = figure();
set(f, 'Units', 'Normalized', 'Position', [0 0 1 1]);

cellid = maprun.mapd03s{get_cell_indices(conerun, offM1)};
urgbs = cell(3,10);
urgbs(3,2:-1:1) = {offurgbs(1) offurgbs(2)};
hist_colors = cell(3,10);
hist_colors{3,2}{1} = colors(1,:);
hist_colors{3,1}{1} = colors(2,:);
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});

cellid = maprun.mapd03s{get_cell_indices(conerun, offM2)};
urgbs = cell(3,10);
urgbs(3,4:-1:3) = {offurgbs(1) offurgbs(2)};
hist_colors = cell(3,10);
hist_colors{3,4}{1} = colors(1,:);
hist_colors{3,3}{1} = colors(2,:);
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});

cellid = maprun.mapd03s{get_cell_indices(conerun, offP)};
urgbs = cell(3,10);
urgbs(3,6:-1:5) = {offurgbs(1) offurgbs(2)};
hist_colors = cell(3,10);
hist_colors{3,6}{1} = colors(1,:);
hist_colors{3,5}{1} = colors(2,:);
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});

cellid = maprun.mapd05s{get_cell_indices(piece.d05s, onM_d05s)};
urgbs = cell(3,10);
urgbs(3,8:-1:7) = {onurgbs(1) onurgbs(2)};
hist_colors = cell(3,10);
hist_colors{3,8}{1} = colors(1,:);
hist_colors{3,7}{1} = colors(2,:);
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});


doubleurgbs = find(urgb.doubles & ~urgb.uds);
intensities = sum(urgb.intensities(cones2plot,doubleurgbs));

cellid = maprun.mapd03s{get_cell_indices(conerun, onP)};
urgbs = cell(3,15);
urgbs(3,14:-1:13) = {onurgbs(1) onurgbs(2)};
urgbs{3,15} = doubleurgbs(intensities == 0.48*6);
hist_colors = cell(3,15);
hist_colors{3,14}{1} = colors(1,:);
hist_colors{3,13}{1} = colors(2,:);
hist_colors{3,15}{1} = [1 0 1];
urgb_raster_subplot(rasterrun, cellid, triggers, urgbs, 'hist_color', hist_colors, rasteropts{:});

sanesubplot(3,5,{1:2 1});
plot_rf_stimmap(conerun, offM1, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);
sanesubplot(3,5,{1:2 2});
plot_rf_stimmap(conerun, offM2, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);
sanesubplot(3,5,{1:2 3});
plot_rf_stimmap(conerun, offP, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);
sanesubplot(3,5,{1:2 4});
plot_rf_stimmap(piece.d05s, onM_d05s, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);
sanesubplot(3,5,{1:2 5});
plot_rf_stimmap(conerun, onP, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, rfopts{:});
axis(bounds);

set(gcf, 'PaperType', 'tabloid');
orient(gcf, 'landscape');

leave(keep_vars{:});
