% Can restore the color/hist_color options to the plotting calls to confirm
% that plot/rasters match up L/R.  These were already setup correctly and
% we decided to go to B/W for this fig.

%% Setup
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
streamedopts = loadopts;
streamedopts.load_sta = true;
streamedopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
streamedopts.set_polarities = {'guess', true};

pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');

keep_vars = {'pieces'; 'loadopts'; 'streamedopts'; 'keep_vars'};


%% Load 2011-12-13-2/data007
piece = '2011-12-13-2';
d01 = load_data([piece '/data001'], streamedopts);
d04s = load_data([piece '/streamed/data004-0/data004-0'], streamedopts);
d04cm = load_data([piece '/data004_data007_data008-norefit/data004-from-data004_data007_data008/data004-from-data004_data007_data008'], streamedopts);
d07cm = load_data([piece '/data004_data007_data008-norefit/data007-from-data004_data007_data008/data007-from-data004_data007_data008'], loadopts);
d04_07_08 = load_data([piece '/data004_data007_data008-norefit/data004_data007_data008-norefit'], loadopts);

d04_07_08.mapd04s = map_ei(d04s, d04_07_08);

d07cm = read_stim_lisp_output(d07cm);
d07cm.stimulus = parse_stim_rgbs(d07cm.stimulus);
d07cm.stimulus = parse_cr_rgbs(d07cm.stimulus);

pieces(piece) = struct('d01', d01, 'd04s', d04s, 'd04cm', d04cm, 'd07cm', d07cm, 'd04_07_08', d04_07_08);
leave(keep_vars{:});


%% Load 2012-09-13-2/data007
piece = '2012-09-13-2';

d00 = load_data([piece '/data000'], streamedopts);
d05s = load_data([piece '/streamed/data005/data005'], streamedopts);
d05 = load_data([piece '/data005'], streamedopts);
d07 = load_data([piece '/data007'], loadopts);
d09 = load_data([piece '/data009'], streamedopts);

d05.mapd05s = map_ei(d05s, d05, 'master_cell_type', 'all');
d07.mapd05s = map_ei(d05s, d07, 'master_cell_type', 'all');
d09.mapd05s = map_ei(d05s, d09, 'master_cell_type', 'all');

d07 = read_stim_lisp_output(d07);
d07.stimulus = parse_stim_rgbs(d07.stimulus);
d07.stimulus = parse_cr_rgbs(d07.stimulus);

pieces(piece) = struct('d00', d00, 'd05s', d05s, 'd05', d05, 'd07', d07, 'd09', d09);
leave(keep_vars{:});

%% Load 2012-09-24-5/data006
piece = '2012-09-24-5';
d00s = load_data([piece '/streamed/data000/data000'], streamedopts);
d03s = load_data([piece '/data003'], streamedopts);
d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);

d03_06_07.mapd03s =  map_ei(d03s, d03_06_07);

d06cm = read_stim_lisp_output(d06cm);
d06cm.stimulus = parse_stim_rgbs(d06cm.stimulus);
d06cm.stimulus = parse_cr_rgbs(d06cm.stimulus);

pieces(piece) = struct('d00s', d00s, 'd03s', d03s, 'd06cm', d06cm, 'd03_06_07', d03_06_07);
leave(keep_vars{:});



%% OFF midgets
% 2011-12-13-2 is nice because has double the number of -0.48 trials (data012 and data013)
% 2012-09-13-2 has very nice cells and more clearly sustained responses
% 2012-09-06-0 has medium RFs, but low baselines and only medium responses
% 2012-09-21-2 has somewhat bigger RFs, but movement is a problem
% 2012-09-24-1 has bigger RFs, but not very clean (data008 better than data006)
% 2012-09-24-5 has somewhat bigger RFs, but offM were not targeted so have to dig around


%% Nice small OFF midget.  Very sustained, very clean RF, 7 cones, highish
% baseline firing (~10 Hz).  Pretty good C/R so might want to use it for that...
% Regions are not quite nicely centered somehow.
piece = pieces('2012-09-13-2');
streamrun = piece.d05s;
streamrun.rgc = 2116;

streamrunnum = get_cell_indices(streamrun, streamrun.rgc);
conerun = piece.d05;
rasterrun = piece.d07;
conerun.rgc = conerun.mapd05s{streamrunnum};
rasterrun.rgc = rasterrun.mapd05s{streamrunnum};
cones2plot = 3:4;
colors = [1 0 0; 0 0 1];

f = figure();

% Find and load in urgbs for raster plot.  It's a 3x2 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = -0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{3,i} = urgbsingles(intensities == intensity);
    hist_colors{3,i}{1} = colors(i,:);
end

% To aid the eye since this is R/L relation in stimmap
urgbs = fliplr(urgbs);
hist_colors = fliplr(hist_colors);

outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', {[0 0 0]}, 'hist_line_width', 2);
sanesubplot(3, 2, {1:2 1:2});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', [0 0 0; 0 0 0], 'fit', false, 'az_pad_factor', 3, 'scaled_up', 10);

set(f, 'position', [1320 204 510 774]);
set(rfax, 'XLim', [88.5 88.5+28], 'YLim', [82.5 82.5+28]);
xlabel(outstructs{3,1}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,2}.AX(2), 'spike rate (Hz)');


% Compare raster to STA timecourse convolved with stimulus
tcrun = piece.d00;
tcrun.rgc = 2116; % Hand matched by EI
tc = tcrun.vision.timecourses(get_cell_indices(tcrun, tcrun.rgc));
period = 1 / 60.35 * 3;
tcx = (1:length(tc.g)) .* period;
stim = zeros(1, length(tc.g));
stim(1:5) = -1;

simtrace  = conv(flipud(tc.g), stim);
simtracex = convx(tcx, tcx);
simtrace = simtrace .* 250;
simtrace = simtrace + 10; % Add baseline
simtrace(simtrace < 0) = 0;

figure;
plot(simtracex, simtrace);


leave(keep_vars{:});


%% Larger off midgets
% 2012-09-24-5: Very nice STAs, responses not as nicely sustained but
% pretty good.  Decent baseline firings.
%   346 Nice, 10 Hz baseline
%   2596 Raster not stable, baseline ~5 Hz. Do not use
%   3677 Nice but baseline firing only ~5 Hz

% Not as clear a sustained response but not bad.  13 strongish cones.
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
conerun.rgc = 346;
rasterrun = piece.d06cm;
maprun = piece.d03_06_07;

conerunnum = get_cell_indices(conerun, conerun.rgc);
rasterrun.rgc = maprun.mapd03s{conerunnum};
cones2plot = 1:2;
colors = [1 0 0; 0 0 1];

f = figure();

% Find and load in urgbs for raster plot.  It's a 3x2 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = -0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{3,i} = urgbsingles(intensities == intensity);
    hist_colors{3,i}{1} = colors(i,:);
end

outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', {[0 0 0]}, 'hist_line_width', 2);
sanesubplot(3, 2, {1:2 1:2});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', [0 0 0; 0 0 0], 'fit', false, 'az_pad_factor', 3, 'scaled_up', 10);
set(f, 'position', [1320 204 510 774]);
set(rfax, 'XLim', [57.5 57.5+28], 'YLim', [160.5 160.5+28]);
set([outstructs{3,1}.AX(2) outstructs{3,2}.AX(2)], 'YLim', [0 70]);
xlabel(outstructs{3,1}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,2}.AX(2), 'spike rate (Hz)');


% Compare raster to STA timecourse convolved with stimulus
tcrun = piece.d00s;
tcrun.rgc = 213; % Hand matched by STA/EI
tc = tcrun.vision.timecourses(get_cell_indices(tcrun, tcrun.rgc));
period = 1 / 60.35 * 3;
tcx = (1:length(tc.g)) .* period;
stim = zeros(1, length(tc.g));
stim(1:5) = -1;

simtrace  = conv(flipud(tc.g), stim);
simtracex = convx(tcx, tcx);
simtrace = simtrace .* 250;
simtrace = simtrace + 10; % Add baseline
simtrace(simtrace < 0) = 0;

figure;
plot(simtracex, simtrace);


leave(keep_vars{:});


%% Nice OFF parasol 1
piece = pieces('2012-09-13-2');
streamrun = piece.d05s;
streamrun.rgc = 5371;

streamrunnum = get_cell_indices(streamrun, streamrun.rgc);
conerun = piece.d05;
rasterrun = piece.d07;
conerun.rgc = conerun.mapd05s{streamrunnum};
rasterrun.rgc = rasterrun.mapd05s{streamrunnum};
cones2plot = 1:2;
colors = [1 0 0; 0 0 1];

f = figure();

% Find and load in urgbs for raster plot.  It's a 3x2 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = -0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{3,i} = urgbsingles(intensities == intensity);
    hist_colors{3,i}{1} = colors(i,:);
end

% To aid the eye since this is R/L relation in stimmap
urgbs = fliplr(urgbs);
hist_colors = fliplr(hist_colors);

outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', {[0 0 0]}, 'hist_line_width', 2);
sanesubplot(3, 2, {1:2 1:2});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', [0 0 0; 0 0 0], 'fit', false, 'az_pad_factor', 3, 'scaled_up', 10);

set(f, 'position', [1320 204 510 774]);
set(rfax, 'XLim', [147.5 147.5+70], 'YLim', [164.5 164.5+70]);
xlabel(outstructs{3,1}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,2}.AX(2), 'spike rate (Hz)');

leave(keep_vars{:});


%% Nice OFF parasol 2
% RF is just a little funky looking.  Another option is 5551, but pretty
% certain we want that one for C/R.  Can try using d09 for RF too, but I
% think actually I like d05 better.
piece = pieces('2012-09-13-2');
streamrun = piece.d05s;
streamrun.rgc = 6841;

streamrunnum = get_cell_indices(streamrun, streamrun.rgc);
conerun = piece.d05;
rasterrun = piece.d07;
conerun.rgc = conerun.mapd05s{streamrunnum};
rasterrun.rgc = rasterrun.mapd05s{streamrunnum};
cones2plot = [2 4];
colors = [1 0 0; 0 0 1];


% Find and load in urgbs for raster plot.  It's a 3x2 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = -0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{3,i} = urgbsingles(intensities == intensity);
    hist_colors{3,i}{1} = colors(i,:);
end

f = figure();
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', {[0 0 0]}, 'hist_line_width', 2);
sanesubplot(3, 2, {1:2 1:2});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', [0 0 0; 0 0 0], 'fit', false, 'az_pad_factor', 3, 'scaled_up', 10);

set(f, 'position', [1320 204 510 774]);
set(rfax, 'XLim', [93.5 93.5+70], 'YLim', [225.5 225.5+70]);
xlabel(outstructs{3,1}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,2}.AX(2), 'spike rate (Hz)');



% Compare raster to STA timecourse convolved with stimulus
tcrun = piece.d00;
tcrun.rgc = 6766; % Hand matched by EI
tc = tcrun.vision.timecourses(get_cell_indices(tcrun, tcrun.rgc));
period = 1 / 60.35 * 3;
tcx = (1:length(tc.g)) .* period;
stim = zeros(1, length(tc.g));
stim(1:5) = -1;

simtrace  = conv(flipud(tc.g), stim);
simtracex = convx(tcx, tcx);
simtrace = simtrace .* 600;
simtrace = simtrace + 7; % Add baseline
simtrace(simtrace < 0) = 0;

figure;
plot(simtracex, simtrace);



leave(keep_vars{:});


%% Nice ON Midget 1
piece = pieces('2011-12-13-2');
streamrun = piece.d04s;
streamrun.rgc = 2821;

streamrunnum = get_cell_indices(streamrun, streamrun.rgc);
conerun = piece.d04cm;
rasterrun = piece.d07cm;
maprun = piece.d04_07_08;
conerun.rgc = maprun.mapd04s{streamrunnum};
rasterrun.rgc = conerun.rgc;
cones2plot = 3:4;
colors = [1 0 0; 0 0 1];

% Find and load in urgbs for raster plot.  It's a 3x2 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = 0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{3,i} = urgbsingles(intensities == intensity);
    hist_colors{3,i}{1} = colors(i,:);
end

% To aid the eye since this is R/L relation in stimmap
urgbs = fliplr(urgbs);
hist_colors = fliplr(hist_colors);

f = figure();
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', {[0 0 0]}, 'hist_line_width', 2);
sanesubplot(3, 2, {1:2 1:2});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', [0 0 0; 0 0 0], 'fit', false, 'az_pad_factor', 3, 'scaled_up', 10);

set(f, 'position', [1320 204 510 774]);
set(rfax, 'XLim', [146.5 146.5+28], 'YLim', [177.5 177.5+28]);
xlabel(outstructs{3,1}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,2}.AX(2), 'spike rate (Hz)');


% Compare raster to STA timecourse convolved with stimulus
tcrun = piece.d01;
tcrun.rgc = 2821; % Hand matched by STA/EI
tc = tcrun.vision.timecourses(get_cell_indices(tcrun, tcrun.rgc));
period = 1 / 60.35 * 1;
tcx = (1:length(tc.g)) .* period;
stim = zeros(1, length(tc.g));
stim(1:15) = 1;

simtrace  = conv(flipud(tc.g), stim);
simtracex = convx(tcx, tcx);
simtrace = simtrace .* 30;
simtrace = simtrace + 1; % Add baseline
simtrace(simtrace < 0) = 0;

figure;
plot(simtracex, simtrace);


leave(keep_vars{:});


%% Nice ON Midget 2
% Regions are not so nice, but otherwise this is the best.  601 is an
% option, but very low baseline firing.
piece = pieces('2011-12-13-2');
streamrun = piece.d04s;
streamrun.rgc = 6601;

streamrunnum = get_cell_indices(streamrun, streamrun.rgc);
conerun = piece.d04cm;
rasterrun = piece.d07cm;
maprun = piece.d04_07_08;
conerun.rgc = maprun.mapd04s{streamrunnum};
rasterrun.rgc = conerun.rgc;
cones2plot = 3:4;
colors = [1 0 0; 0 0 1];

% Find and load in urgbs for raster plot.  It's a 3x2 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = 0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{3,i} = urgbsingles(intensities == intensity);
    hist_colors{3,i}{1} = colors(i,:);
end

% To aid the eye since this is R/L relation in stimmap
urgbs = fliplr(urgbs);
hist_colors = fliplr(hist_colors);


f = figure();
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'hist_color', hist_colors, 'hist_line_width', 2);
sanesubplot(3, 2, {1:2 1:2});
plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, 'fit', false, 'az_pad_factor', 3, 'scaled_up', 10);

set(f, 'position', [1320 204 510 774]);
xlabel(outstructs{3,1}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,2}.AX(2), 'spike rate (Hz)');

leave(keep_vars{:});


%% ON Parasol
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
conerun.rgc = 7237;
rasterrun = piece.d06cm;
maprun = piece.d03_06_07;

conerunnum = get_cell_indices(conerun, conerun.rgc);
rasterrun.rgc = maprun.mapd03s{conerunnum};
cones2plot = [2 4];
colors = [1 0 0; 0 0 1];


% Find and load in urgbs for raster plot.  It's a 3x2 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = 0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{3,i} = urgbsingles(intensities == intensity);
    hist_colors{3,i}{1} = colors(i,:);
end

% To aid the eye since this is R/L relation in stimmap
urgbs = fliplr(urgbs);
hist_colors = fliplr(hist_colors);

% Add double
urgbdoubles = find(urgb.doubles & ~urgb.uds);
intensities = sum(urgb.intensities(cones2plot, urgbdoubles));
urgbs{3,3} = urgbdoubles(intensities == intensity*2);
hist_colors{3,3}{1} = [1 0 1];


f = figure();
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'hist_color', hist_colors, 'hist_line_width', 2);
sanesubplot(3, 3, {1:2 1:2});
plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, 'fit', false, 'az_pad_factor', 4.5, 'scaled_up', 10);

set(f, 'position', [1320 204 510 774]);
xlabel(outstructs{3,1}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,3}.AX(2), 'spike rate (Hz)');

leave(keep_vars{:});


%% Searching...
d = pieces('2012-09-24-5');

% Pick trials (same d04cm d06cm)
urgb = d.d06cm.stimulus.urgb;
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
ploturgbs = [num2cell(fliplr(singleon)); num2cell(singleoff)];

% Colors to match regions
colors = jet(4);
colors = mat2cell(colors, [1 1 1 1], 3);
colors = num2cell(colors);
colors = [colors'; colors'];

conerun = d.d03s;
rasterrun = d.d06cm;
map = d.d03_06_07.mapd03s;
for cellid = conerun.cell_types{5}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'add_width', 2, 'start', -0.2, 'color', {[0.5 0.5 0.5]}, 'hist_color', colors, 'hist_line_width', 2);

    sanesubplot(1, 3, {1 3});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end

leave(keep_vars{:});


%% Searching...
d = pieces('2012-09-13-2');

% Pick trials (same d04cm d06cm)
urgb = d.d07.stimulus.urgb;
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
ploturgbs = [num2cell(fliplr(singleon)); num2cell(singleoff)];

% Colors to match regions
colors = jet(4);
colors = mat2cell(colors, [1 1 1 1], 3);
colors = num2cell(colors);
colors = [colors'; colors'];

conerun = d.d05s;
rasterrun = d.d07;
map = d.d07.mapd05s;
for cellid = conerun.cell_types{10}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'add_width', 2, 'start', -0.2, 'color', {[0.5 0.5 0.5]}, 'hist_color', colors, 'hist_line_width', 2);

    sanesubplot(1, 3, {1 3});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end

leave(keep_vars{:});
