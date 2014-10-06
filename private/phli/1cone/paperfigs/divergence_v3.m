%% Setup
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'pieces'; 'loadopts'; 'staopts'; 'keep_vars'};


%% Custom plot to get best rasters and RFs
% 2012-09-24-5 M-, M+, P- & surround M-
piece = '2012-09-24-5';
d01s = load_data([piece '/data001'], staopts);
d03s = load_data([piece '/data003'], staopts);
d04s = load_data([piece '/data004'], loadopts);
d05s = load_data([piece '/data005'], staopts);
d03cm = load_data([piece '/d03-06-07-norefit/data003/data003'], staopts);
d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], staopts);
d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d04cm = load_data([piece '/d01-04-05-norefit/data004/data004'], loadopts);
d01_04_05 = load_data([piece '/d01-04-05-norefit/d01-04-05-norefit'], loadopts);
d04s.mapd01s = map_ei(d01s, d04s);
d04s.mapd03s = map_ei(d03s, d04s);
d03_06_07.mapd03s = map_ei(d03s, d03_06_07);
d03_06_07.mapd05s = map_ei(d05s, d03_06_07);
d01_04_05.mapd03s = map_ei(d03s, d01_04_05);
d04s = read_stim_lisp_output(d04s);
d04s.stimulus = parse_stim_rgbs(d04s.stimulus);
d04s.stimulus = parse_cr_rgbs(d04s.stimulus);
d06cm = read_stim_lisp_output(d06cm);
d06cm.stimulus = parse_stim_rgbs(d06cm.stimulus);
d06cm.stimulus = parse_cr_rgbs(d06cm.stimulus);

%%
cone2plot = 4;

%%
x = 107.5; y = 70.5; del = 50; 
bounds = [x x+del y y+del];

% Override complex sigstix in favor of rstd threshold
thresh = 3.5;
datarun = d01s;
rgc = 2852;
sta = get_sta(datarun, rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
datarun.stas.rfs{datarun.cell_nums(rgc)} = [];
datarun.stas.marks{datarun.cell_nums(rgc)} = simplemarks;
sanesubplot(3, 8, {1:2 1:2});
plot_rf_stimmap(datarun, rgc, d04s.stimulus.mapnycpoly{1}(cone2plot), 'fit', false, 'colors', 'k', 'scaled_up', 10, 'bounds', bounds);

thresh = 3.5;
datarun = d03s;
rgc = 2855;
sta = get_sta(datarun, rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
datarun.stas.rfs{datarun.cell_nums(rgc)} = [];
datarun.stas.marks{datarun.cell_nums(rgc)} = simplemarks;
sanesubplot(3, 8, {1:2 3:4});
plot_rf_stimmap(datarun, rgc, d04s.stimulus.mapnycpoly{1}(cone2plot), 'fit', false, 'colors', 'k', 'scaled_up', 10, 'bounds', bounds);

thresh = 3.5; % Marks seem better at 5, but RF looks as good here
datarun = d01s;
rgc = 2521;
sta = get_sta(datarun, rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
datarun.stas.rfs{datarun.cell_nums(rgc)} = [];
datarun.stas.marks{datarun.cell_nums(rgc)} = simplemarks;
sanesubplot(3, 8, {1:2 5:6});
plot_rf_stimmap(datarun, rgc, d04s.stimulus.mapnycpoly{1}(cone2plot), 'fit', false, 'colors', 'k', 'scaled_up', 10, 'bounds', bounds);

thresh = 3.5;
datarun = d01s;
rgc = 2642;
sta = get_sta(datarun, rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
datarun.stas.rfs{datarun.cell_nums(rgc)} = [];
datarun.stas.marks{datarun.cell_nums(rgc)} = simplemarks;
sanesubplot(3, 8, {1:2 7:8});
plot_rf_stimmap(datarun, rgc, d04s.stimulus.mapnycpoly{1}(cone2plot), 'fit', false, 'colors', 'k', 'scaled_up', 10, 'bounds', bounds);


%%
start = -0.1;
stop = 0.75;

d04s.urgb = d04s.stimulus.urgb;
d04s.offi = d04s.urgb.singles & d04s.urgb.intensities(cone2plot,:) == -1.44;
d04s.oni = d04s.urgb.singles & d04s.urgb.intensities(cone2plot,:) == 1.44;
d04cm.stimulus = d04s.stimulus;
d04cm.urgb = d04s.stimulus.urgb;
d04cm.offi = d04s.offi;
d04cm.oni = d04s.oni;
d06cm.urgb = d06cm.stimulus.urgb;
d06cm.offi = d06cm.urgb.singles & d06cm.urgb.intensities(cone2plot,:) == -1.44;
d06cm.oni = d06cm.urgb.singles & d06cm.urgb.intensities(cone2plot,:) == 1.44;

urgbs = cell(3, 8);
urgbs{3,1} = {find(d04s.offi)};
outstruct = urgb_raster_subplot(d04s, d04s.mapd01s{get_cell_indices(d01s, 2852)}, d04s.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k', 'hist_line_width', 2, 'start', start, 'stop', stop);
ax = outstruct{3,1}.AX;
axis(ax(1), [-0.08 0.75 0 41]);
axis(ax(2), [-0.08 0.75 0 60]);
pbaspect(ax(1), [4 3 1]);
pbaspect(ax(2), [4 3 1]);

urgbs = cell(3, 8);
urgbs{3,3} = {find(d04cm.oni)};
outstruct = urgb_raster_subplot(d04cm, d01_04_05.mapd03s{get_cell_indices(d03s, 2855)}, d04cm.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k', 'hist_line_width', 2, 'start', start, 'stop', stop);
ax = outstruct{3,3}.AX;
axis(ax(1), [-0.08 0.75 0 41]);
axis(ax(2), [-0.08 0.75 0 15]);
pbaspect(ax(1), [4 3 1]);
pbaspect(ax(2), [4 3 1]);

% urgbs = cell(3, 8);
% urgbs{3,3} = {find(d06cm.oni)};
% outstruct = urgb_raster_subplot(d06cm, d03_06_07.mapd03s{get_cell_indices(d03s, 2855)}, d06cm.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k', 'hist_line_width', 2, 'start', start, 'stop', stop);
% ax = outstruct{3,3}.AX;
% axis(ax(1), [-0.08 0.75 0 41]);
% axis(ax(2), [-0.08 0.75 0 15]);
% pbaspect(ax(1), [4 3 1]);
% pbaspect(ax(2), [4 3 1]);

urgbs = cell(3, 8);
urgbs{3,5} = {find(d04s.offi)};
outstruct = urgb_raster_subplot(d04s, d04s.mapd01s{get_cell_indices(d01s, 2521)}, d04s.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k', 'hist_line_width', 2, 'start', start, 'stop', stop);
ax = outstruct{3,5}.AX;
axis(ax(1), [-0.08 0.75 0 41]);
axis(ax(2), [-0.08 0.75 0 40]);
pbaspect(ax(1), [4 3 1]);
pbaspect(ax(2), [4 3 1]);

urgbs = cell(3, 8);
urgbs{3,7} = {find(d04s.offi)};
outstruct = urgb_raster_subplot(d04s, d04s.mapd01s{get_cell_indices(d01s, 2642)}, d04s.triggers(1:2:end), urgbs, 'color', {[0.5 0.5 0.5]}, 'hist_color', 'k', 'hist_line_width', 2, 'start', start, 'stop', stop);
ax = outstruct{3,7}.AX;
axis(ax(1), [-0.08 0.75 0 41]);
axis(ax(2), [-0.08 0.75 0 20]);
pbaspect(ax(1), [4 3 1]);
pbaspect(ax(2), [4 3 1]);


%% Surround versus baseline discriminability
start = 0.05;
stop = 0.25;

res1 = outstruct{3,3}.res;
res1 = res1(res1(:,1) > start & res1(:,1) < stop,:);
spikecount1 = zeros(max(res1(:,2)),1);
for i = 1:size(res1,1)
    trialnum = res1(i,2);
    spikecount1(trialnum) = spikecount1(trialnum)+1;
end

spikecount0 = zeros(max(res0(:,2)),1);

mu1 = mean(spikecount1);
mu2 = mean(spikecount2);
mu0 = mean(spikecount0);
sigma1 = std(spikecount1);
sigma2 = std(spikecount2);
sigma0 = std(spikecount0);

d1 = (mu1 - mu0) / sigma0
d2 = (mu2 - mu0) / sigma0



%% Old stuff below


%% Load 2012-09-24-5
pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');
piece = '2012-09-24-5';
d01s = load_data([piece '/data001'], staopts);
d03s = load_data([piece '/data003'], staopts);
d04s = load_data([piece '/data004'], loadopts);
d05s = load_data([piece '/data005'], staopts);
d03cm = load_data([piece '/d03-06-07-norefit/data003/data003'], staopts);
d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], staopts);
d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d04s.mapd01s = map_ei(d01s, d04s);
d04s.mapd03s = map_ei(d03s, d04s);
d03_06_07.mapd03s = map_ei(d03s, d03_06_07);
d03_06_07.mapd05s = map_ei(d05s, d03_06_07);
d04s = read_stim_lisp_output(d04s);
d04s.stimulus = parse_stim_rgbs(d04s.stimulus);
d04s.stimulus = parse_cr_rgbs(d04s.stimulus);
d06cm = read_stim_lisp_output(d06cm);
d06cm.stimulus = parse_stim_rgbs(d06cm.stimulus);
d06cm.stimulus = parse_cr_rgbs(d06cm.stimulus);
pieces(piece) = struct('d01s', d01s, 'd03s', d03s, 'd04s', d04s, 'd05s', d05s, 'd03cm', d03cm, 'd06cm', d06cm, 'd07cm', d07cm, 'd03_06_07', d03_06_07);
leave(keep_vars{:});


%% Divergence 1, 2012-09-24-5, Ms & OFFP
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d06cm;
map = piece.d03_06_07.mapd03s;

% Possible surround OFF midgets
%surrounds = [1336 2176 2251 2252 2341 2401 2431 2492 2507 2566 2641 2662 2686 2776 2851 2866 3121 3151 3152 3196 3226 3316 3452 3456 4216];
%rgcs = surrounds(21:25);
% Great: 2776
% Nice edge: 2641, 3152
% Good farther: 2866, 3151
% Good: 2492
% Okay: 2401
% Okay farther: 3121, 3226

rgcs = [];
rgcs(end+1) = 2851; % offM
rgcs(end+1) = 2855; % onM
rgcs(end+1) = 2521; % offP
rgcs(end+1) = 2776; % surround OFFM

cones2plot = 4;
bounds = [105 160 70 125];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax(:,[1 3 5]), 'YTick', []);

cell2mat(cellfun(@(S)(S.AX'), outstructs, 'UniformOutput', false))


%% Divergence 1.1, 2012-09-24-5, Ms & OFFP (rasterrun with fewer regions)
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d04s;
map = rasterrun.mapd03s;

rgcs = [];
rgcs(end+1) = 2851; % offM
rgcs(end+1) = 2855; % onM
rgcs(end+1) = 2521; % offP
rgcs(end+1) = 2776; % surround OFFM

cones2plot = 4;
bounds = [105 160 70 125];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax(:,[1 3 5]), 'YTick', []);

ax = cell2mat(cellfun(@(S)(S.AX'), outstructs, 'UniformOutput', false));


%% Divergence 1.2, 2012-09-24-5, Ms & OFFP (rasterrun with fewer regions, different surround)
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d04s;
map = rasterrun.mapd03s;

rgcs = [];
rgcs(end+1) = 2851; % offM           % Check d01 id2852
rgcs(end+1) = 2855; % onM
rgcs(end+1) = 2521; % offP           % Check d01 id 2521
rgcs(end+1) = 2507; % surround OFFM  % Looks better in d01 id2642

cones2plot = 4;
bounds = [105 160 70 125];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax(:,[1 3 5]), 'YTick', []);

ax = cell2mat(cellfun(@(S)(S.AX'), outstructs, 'UniformOutput', false));


%% Cleanup
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax(:,2:end), 'YTick', []);
set(rax(:,1:end-1), 'YTick', []);

set(gcf, 'PaperType', 'tabloid');
orient(gcf, 'landscape');

% Match Y lims across whole fig?
% setaxesy(cellfun(@(s)(s.AX(2)), outstructs))