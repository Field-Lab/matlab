%% Look at covariance of divergence cells
% Possible figure for 1cone/elementary paper.  But results not looking
% strong so probably drop it.

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

%% Load 2012-09-13-2
piece = '2012-09-13-2';
d01s = load_data([piece '/streamed/data001/data001'], staopts);
d05s = load_data([piece '/streamed/data005/data005'], staopts);
d04cm = load_data([piece '/d01-04-05-norefit/data004/data004'], loadopts);
d07 = load_data([piece '/data007'], loadopts);
d04cm.mapd01s = map_ei(d01s, d04cm);
d07.mapd05s = map_ei(d05s, d07);
d04cm = read_stim_lisp_output(d04cm);
d04cm.stimulus = parse_stim_rgbs(d04cm.stimulus);
d04cm.stimulus = parse_cr_rgbs(d04cm.stimulus);
d07 = read_stim_lisp_output(d07);
d07.stimulus = parse_stim_rgbs(d07.stimulus);
d07.stimulus = parse_cr_rgbs(d07.stimulus);
pieces(piece) = struct('d01s', d01s, 'd05s', d05s, 'd04cm', d04cm, 'd07', d07);
leave(keep_vars{:});

%% Kludge: plot, but actually just interested in rasters so close plot
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d06cm;
map = piece.d03_06_07.mapd03s;
rgcs = [];
rgcs(end+1) = 2851; % offM
rgcs(end+1) = 2521; % offP
cones2plot = 4;
bounds = [110 160 70 120];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
close(gcf);

piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d06cm;
maprun = piece.d03_06_07;
rgcs = [];
rgcs(end+1) = 3646; % offM
rgcs(end+1) = 3616; % offP
cones2plot = 3;
bounds = [160 220 15 85];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
close(gcf);

piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d06cm;
map = piece.d03_06_07.mapd03s;
rgcs = [];
rgcs(end+1) = 3677; % offM
rgcs(end+1) = 4172; % offP
cones2plot = 2;
bounds = [190 240 75 140];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
close(gcf);

piece = pieces('2012-09-13-2');
conerun = piece.d01s;
rasterrun = piece.d04cm;
map = piece.d04cm.mapd01s;
rgcs = [];
rgcs(end+1) = 7516; % offM
rgcs(end+1) = 7666; % offP
cones2plot = 1;
bounds = [80 150 230 290];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);

piece = pieces('2012-09-13-2');
conerun = piece.d01s;
rasterrun = piece.d04cm;
map = piece.d04cm.mapd01s;
rgcs = [];
rgcs(end+1) = 7516; % offM
rgcs(end+1) = 6993; % offP
cones2plot = 1;
bounds = [80 150 230 290];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);

piece = pieces('2012-09-13-2');
conerun = piece.d01s;
rasterrun = piece.d04cm;
map = piece.d04cm.mapd01s;
rgcs = [];
rgcs(end+1) = 6993; % offP
rgcs(end+1) = 7666; % offP
cones2plot = 1;
bounds = [80 150 230 290];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);

piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 2386; % offM
rgcs(end+1) = 2506; % offP
cones2plot = 2;
bounds = [135 190 65 120];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);

piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 2386; % offM
rgcs(end+1) = 2506; % offP
cones2plot = 4;
bounds = [135 190 65 120];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);

piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 6991; % offM
rgcs(end+1) = 6601; % offP
cones2plot = 1;
bounds = [105 185 230 285];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);

piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 6991; % offM
rgcs(end+1) = 6841; % offP
cones2plot = 1;
bounds = [105 185 230 285];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);

piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 6841; % offP
rgcs(end+1) = 6601; % offP
cones2plot = 1;
bounds = [105 185 230 285];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);

%% Get params for each trial for OFF M and OFF P
res1 = outstructs{2}.res;
res2 = outstructs{4}.res;
ntrials = 40;

res = res1;
ttfs = NaN(1,ntrials);
sums = zeros(1,ntrials);
for i = 1:ntrials
    spikes = res(res(:,2) == i,1);
    spikes = spikes(spikes > 0);
    if isempty(spikes), continue; end

    ttfs(i) = min(spikes);
    sums(i) = length(spikes);
end
ttfs1 = ttfs;
sums1 = sums;

res = res2;
ttfs = NaN(1,ntrials);
sums = zeros(1,ntrials);
for i = 1:ntrials
    spikes = res(res(:,2) == i,1);
    spikes = spikes(spikes > 0);
    if isempty(spikes), continue; end

    ttfs(i) = min(spikes);
    sums(i) = length(spikes);
end
ttfs2 = ttfs;
sums2 = sums;

realcorr = corr(ttfs1', ttfs2');

figure
plot(ttfs1,ttfs2,'.');
axis equal
title(sprintf('r = %f', realcorr));

%% Permutation test

allpermcorrs = [];

% Repeat this part as desired.
nperms = 10000;
permcorrs = zeros(1,nperms);
for i = 1:nperms;
    permed = ttfs2(randperm(ntrials));
    permcorrs(i) = corr(ttfs1', permed');
end
allpermcorrs = [allpermcorrs permcorrs];

p = sum(allpermcorrs > realcorr) / length(allpermcorrs);

%% Plot permutation test results
figure
hist(allpermcorrs, 100);
hold on;
maxY = max(get(gca, 'YLim'));
[arrX arrY] = dsxy2figxy([realcorr realcorr], [maxY/2 maxY/50]);
annotation('arrow', arrX, arrY);
title(sprintf('p = %f', p));