% Basics
piece = '2012-09-24-5';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
streamedopts = loadopts;
streamedopts.load_sta = true;
streamedopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
streamedopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'streamedopts'; 'keep_vars'; 'd'};


%% Load data
d.d01s = load_data([piece '/data001'], streamedopts);
d.d01_04_05 = load_data([piece '/d01-04-05-norefit/d01-04-05-norefit'], loadopts);
d.d04cm = load_data([piece '/d01-04-05-norefit/data004/data004'], loadopts);
d.d05cm = load_data([piece '/d01-04-05-norefit/data005/data005'], loadopts);
d.d01_04_05.mapd01s = map_ei(d.d01s, d.d01_04_05);
d.d04cm = read_stim_lisp_output(d.d04cm);
d.d04cm.stimulus = parse_stim_rgbs(d.d04cm.stimulus);

d.d03s = load_data([piece '/data003'], streamedopts);
d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts);
d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);
d.d03_06_07.mapd03s = map_ei(d.d03s, d.d03_06_07);


%% Surround tangential
conerun = d.d03s;
conerun.rgc = 2776;
rasterrun = d.d06cm;
map = d.d03_06_07.mapd03s;

conerunnum = get_cell_indices(conerun, conerun.rgc);
rasterrun.rgc = map{conerunnum};
cones2plot = [4 3 1];
colors = jet(length(cones2plot));

% Find and load in urgbs for raster plot.  It's a Nx3 cell array with the
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
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'start', -0.1, 'stop', 0.7, 'color', {[0.5 0.5 0.5]}, 'hist_color', {[0 0 0]}, 'hist_line_width', 2);
sanesubplot(3, 1, {1:2 1});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', zeros(4,3), 'fit', false, 'az_pad_factor', 10, 'scaled_up', 10);

set(f, 'position', [1965           4        1405         990]);
set(rfax, 'XLim', [127 127+40], 'YLim', [84 84+40]);
xlabel(outstructs{3,2}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,3}.AX(2), 'spike rate (Hz)');

leave(keep_vars{:});


%% Series, 1: Edge
conerun = d.d01s;
conerun.rgc = 871;
rasterrun = d.d04cm;
map = d.d01_04_05.mapd01s;

conerunnum = get_cell_indices(conerun, conerun.rgc);
rasterrun.rgc = map{conerunnum};
cones2plot = [4 3 1];
colors = jet(length(cones2plot));

% Find and load in urgbs for raster plot.  It's a Nx3 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = -0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{i,3} = urgbsingles(intensities == intensity);
    hist_colors{i,3}{1} = colors(i,:);
end

f = figure();
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', {[0 0 0]}, 'hist_line_width', 2);
sanesubplot(3, 3, {1:3 1:2});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', [0 0 0; 0 0 0; 0 0 0], 'fit', false, 'az_pad_factor', 10, 'scaled_up', 10);

set(f, 'position', [1311 194 1283 774]);
set(rfax, 'Xlim', [85 85+45], 'Ylim', [125 125+45]);
set([outstructs{2,3}.AX(2) outstructs{3,3}.AX(2)], 'YLim', [0 30]);
xlabel(outstructs{3,3}.AX(1), 'time (s)');
ylabel(outstructs{2,3}.AX(1), 'trial #');
ylabel(outstructs{2,3}.AX(2), 'spike rate (Hz)');

leave(keep_vars{:});

%% Series, 2: Surround
conerun = d.d01s;
conerun.rgc = 796;
rasterrun = d.d04cm;
map = d.d01_04_05.mapd01s;

conerunnum = get_cell_indices(conerun, conerun.rgc);
rasterrun.rgc = map{conerunnum};
cones2plot = [4 3 1];
colors = jet(length(cones2plot));

% Find and load in urgbs for raster plot.  It's a Nx3 cell array with the
% top 2 rows blank (leaving room for the RF stimmap).
urgb = rasterrun.stimulus.urgb;
urgbsingles = find(urgb.singles);
intensity = -0.48*3;
for i = 1:length(cones2plot)
    cone = cones2plot(i);
    intensities = urgb.intensities(cone, urgb.singles);
    urgbs{i,3} = urgbsingles(intensities == intensity);
    hist_colors{i,3}{1} = colors(i,:);
end

f = figure();
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', {[0 0 0]}, 'hist_line_width', 2);
sanesubplot(3, 3, {1:3 1:2});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', [0 0 0; 0 0 0; 0 0 0], 'fit', false, 'az_pad_factor', 10, 'scaled_up', 10);

set(f, 'position', [1311 194 1283 774]);
set(rfax, 'Xlim', [85 85+45], 'Ylim', [125 125+45]);
xlabel(outstructs{3,3}.AX(1), 'time (s)');
ylabel(outstructs{2,3}.AX(1), 'trial #');
ylabel(outstructs{2,3}.AX(2), 'spike rate (Hz)');

leave(keep_vars{:});

%% Surround single
conerun = d.d01s;
conerun.rgc = 3316;
rasterrun = d.d04cm;
map = d.d01_04_05.mapd01s;

conerunnum = get_cell_indices(conerun, conerun.rgc);
rasterrun.rgc = map{conerunnum};
cones2plot = [4 3 2];
colors = jet(length(cones2plot));

% Find and load in urgbs for raster plot.  It's a Nx3 cell array with the
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
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, rasterrun.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', hist_colors, 'hist_line_width', 2);
sanesubplot(3, 1, {1:2 1});
plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', colors, 'fit', false, 'az_pad_factor', 10, 'scaled_up', 10);

set(f, 'position', [ 1974           3         658         983]);
xlabel(outstructs{3,2}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,3}.AX(2), 'spike rate (Hz)');

leave(keep_vars{:});

%%  Cusp single
conerun = d.d03s;
conerun.rgc = 5071;
rasterrun = d.d06cm;
map = d.d03_06_07.mapd03s;

conerunnum = get_cell_indices(conerun, conerun.rgc);
rasterrun.rgc = map{conerunnum};
cones2plot = [4 1 2 3];
colors = jet(length(cones2plot));

% Find and load in urgbs for raster plot.  It's a Nx3 cell array with the
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
sanesubplot(3, 1, {1:2 1});
rfax = plot_rf_stimmap(conerun, conerun.rgc, rasterrun.stimulus.mapnycpoly{1}(cones2plot), 'colors', zeros(4,3), 'fit', false, 'az_pad_factor', 10, 'scaled_up', 10);

set(f, 'position', [1965           4        1405         990]);
set(rfax, 'XLim', [212 212+35], 'YLim', [140 140+35]);
xlabel(outstructs{3,2}.AX(1), 'time (s)');
ylabel(outstructs{3,1}.AX(1), 'trial #');
ylabel(outstructs{3,3}.AX(2), 'spike rate (Hz)');

leave(keep_vars{:});

%% Try to find surround in RF; spatial profile 
conerun = d.d03s;
conerun.rgc = 5071;
conerun = load_cones(conerun);
coneweightopts = struct('thresh', 0, 'radius', [0 8], 'polarity', -1, 'contiguity', false,'scale', 3.0);
[mosaic_weights, selection, extras] = select_cone_weights(conerun, conerun.rgc, coneweightopts);

% Combine STAs for better SNR
sta1 = get_sta(d.d01s, 5071);
sta3 = get_sta(d.d03s, 5071);
sta = (sta1 + sta3) ./ 2;
rf = rf_from_sta(sta);

% Upscale and average over surround cone centers
roi = rf(75:225,150:300);
roiscale = 10;
roi = matrix_scaled_up(roi, roiscale);

% Get new cone centers
ncones = size(conerun.cones.centers,1);
scaledconecenters = ((conerun.cones.centers - repmat([149 74], [ncones 1])) - 0.5) * roiscale + 0.5;

% Check
imshow(norm_image(-roi));
hold on
plot(scaledconecenters(selection,1), scaledconecenters(selection,2), '.y');

% Calculate spatial average
r = 200;
avsurrcone = zeros(r*2 + 1);
for i = find(selection)'
    ctr = round(scaledconecenters(i,:));
    x = (-r:r) + ctr(1);
    y = (-r:r) + ctr(2);
    currroi = roi(x,y);
%    figure; imagesc(currroi);
    avsurrcone = avsurrcone + roi(x,y);
end
figure; imagesc(-avsurrcone); colormap gray

%% Try to find surround in RF
conerun = d.d03s;
conerun.rgc = 5071;
conerun = load_cones(conerun);
coneweightopts = struct('thresh', 0, 'radius', [0 10], 'polarity', 0, 'contiguity', false,'scale', 3.0);
[mosaic_weights, selection, extras] = select_cone_weights(conerun, conerun.rgc, coneweightopts);

cellnum = get_cell_indices(conerun, conerun.rgc);
rgcfit = conerun.stas.fits{cellnum};
selectedconecenters = conerun.cones.centers(selection,:);
del = selectedconecenters - repmat(rgcfit.mean, [sum(selection) 1]);
del2 = del.^2;
distances = sqrt(sum(del2, 2));

plot(distances, mosaic_weights(selection), '.');
hold on;
plot ([0 max(distances)], [0 0], '--');

%% Searching
%% Pick trials (same d04cm d06cm)
urgb = d.d04cm.stimulus.urgb;
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
ploturgbs = [num2cell(fliplr(singleon)); num2cell(singleoff)];

%% Colors to match regions
colors = jet(4);
colors = mat2cell(colors, [1 1 1 1], 3);
colors = num2cell(colors);
colors = [colors'; colors'];

%% Picks d04cm
conerun = d.d01s;
rasterrun = d.d04cm;
map = d.d01_04_05.mapd01s;

% Edge cells: [332 335 871 2192 2431 2584 2641 2642 3676 5626 6256 6391 6392 6496 6511 6601 6631 6664 6841]
% Surround cells: [122 138 347 423 529 707 828 2747 2776 3056 3736 3797 3811 3828 4051 4112 4127 4173 6466 7711]
% Surround cells: [451 796 886 916 1051 1081 1111 1338 1727 1741 2356 3316 3632 3646 3721 4142 4143 4338 4411 4771 5209 5237 5492 5506 5582 5716 5986 6272 6560 6707 6811 6856 6991 7081 7096 7247 7248 7264 7366 7396 7561 7578]
cuspeffect = [6391 6496 6392    6560];
distanceeffect = [5626 6256     871 2192    1727 5237 796 916 1051 7366 7651];
justedge = [451 2431];
surround = [3316    529 347];

%% Picks d06cm
conerun = d.d03s;
rasterrun = d.d06cm;
map = d.d03_06_07.mapd03s;

% Edge cells: 77 421 466 856 1051 1111 2641 2851 3152 3211 4036 4141 4336 5071 5191
% Surround:  48 122 138 961 1217 1291 1382 2401 3811 4066 4186 4742 5266 5326 5522 5536 5732 6706 7081 7336 7351 7547 7561 7576 7681 7698
edge = [5071 1051];
distanceeffect = [421 7698 7561 6706];
cuspeffect = 5326;

%% Plotting
for cellid = [3031, 2866, 2776, 2686, 2492, 1336, 1186, 7081, 7351, 7576, 7561, 6706, 3031, 2492, 1186, 452, 451, 332, 5732, 5536, 5522, 5326, 5266, 5056, 4846, 4742, 4576, 3811, 3826, 3721, 3676, 3226, 7698, 7681, 7547 7576, 7561, 138, 122, 48, 827, 556, 541, 1111, 1051, 466, 751, 7576, 1382, 1291, 1217, 961, 1741, 1546, 1336, 1276, 1186, 827, 556, 541, 452, 451, 332, 1488, 1353, 1231, 1081, 3226, 4742, 4186, 5056, 4846, 4742, 4576, 4516 4339, 4337, 4306, 4051, 4066, 4576, 4516, 4411, 4339, 4337, 4306, 4051, 3826, 3721, 3691, 3676];
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'add_width', 2, 'start', -0.2, 'color', {[0.5 0.5 0.5]}, 'hist_color', colors, 'hist_line_width', 2);

    sanesubplot(1, 3, {1 3});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end