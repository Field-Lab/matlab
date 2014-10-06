%% Setup
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');
keep_vars = {'pieces'; 'loadopts'; 'staopts'; 'keep_vars'};

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

%% Load 2012-09-24-5
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




%% Divergence 1b, 2012-09-24-5, many M-
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d04s;
map = rasterrun.mapd03s;

% Possible surround OFF midgets
% surrounds = [1336 2176 2251 2252 2341 2401 2431 2492 2507 2566 2641 2662 2686 2776 2851 2866 3121 3151 3152 3196 3226 3316 3452 3456 4216];
% rgcs = surrounds(21:25);

rgcs = [];
rgcs = [rgcs 2851 3152]; % center
rgcs = [rgcs 2507 2641]; % edge
rgcs = [rgcs 2251 2252 2341 2401 2566 2686 2776 2866 3151 3196 3316 3452 4216]; %surround

cones2plot = 4;
bounds = [90 195 45 150];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs(1:8), cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);
set(lax(:,2:end),   'YTick', []);
set(rax(:,[1 3 5]), 'YTick', []);

cones2plot = 4;
bounds = [95 200 45 150];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs(9:end), cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);
set(lax(:,2:end),   'YTick', []);
set(rax(:,[1 3 5]), 'YTick', []);


%% Divergence 1c
piece = pieces('2012-09-24-5');
conerun = piece.d01s;
rasterrun = piece.d04s;
map = rasterrun.mapd01s;

% surrounds = [122 138 261 332 335 336 451 466 529 796 856 871 886 916 1051 1111 1173 1186 1188 1324 1471 2492 7187 7232 7366 7561 7578];
% rgcs = surrounds(21:27);


rgcs = [];
rgcs = [rgcs 335 466]; % center M-
rgcs = [rgcs 332 451 871]; % edge
rgcs = [rgcs 122 529 796 856 886 916 1051 1186 1471 2492 7366 7561 7578]; % surround
% rgcs = [rgcs 138 1173 7187]; % very weak surround
% rgcs = [rgcs 261 1188]; % center P-
% rgcs = [rgcs 1324]; % surround P-

cones2plot = [4 1];
bounds = [45 155 100 210];

[rfax outstructs] = divergence_compound_plot(conerun, rasterrun, map, rgcs(1:6), cones2plot, bounds, 'cmf', @winter, 'rasteropts', {'setaxesy', false});
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);
set(lax(:,2:end),   'YTick', []);
set(rax(:,[1 3 5]), 'YTick', []);

[rfax outstructs] = divergence_compound_plot(conerun, rasterrun, map, rgcs(7:12), cones2plot, bounds, 'cmf', @winter, 'rasteropts', {'setaxesy', false});
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);
set(lax(:,2:end),   'YTick', []);
set(rax(:,[1 3 5]), 'YTick', []);

[rfax outstructs] = divergence_compound_plot(conerun, rasterrun, map, rgcs(13:end), cones2plot, bounds, 'cmf', @winter, 'rasteropts', {'setaxesy', false});
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);
set(lax(:,2:end),   'YTick', []);
set(rax(:,[1 3 5]), 'YTick', []);


%% Divergence 1d
piece = pieces('2012-09-24-5');
conerun = piece.d01s;
rasterrun = piece.d04s;
map = rasterrun.mapd01s;

surrounds = [1202 2822 3031 3046 4411 4486 4501 4771 4966 5237 5492 5506 5509 5626 5716 6062 6256 6272 6286 6287 6302 6391 6392 6423 6558 6560 6587 6707 6841 6856 6858 6886 6917 6931 7081 7173];



%% Divergence 2, 2012-09-24-5 OFF M and P
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
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax(:,1:2:end), 'YTick', []);


%% Divergence 3, 2012-09-24-5, ON & OFF M
piece = pieces('2012-09-24-5');
conerun = piece.d03s;
rasterrun = piece.d06cm;
map = piece.d03_06_07.mapd03s;
rgcs = [];
rgcs(end+1) = 4591; % offM
rgcs(end+1) = 4597; % onM
cones2plot = 3;
bounds = [235 270 105 140];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax(:,1:2:end), 'YTick', []);


%% Divergence 4, 2012-09-24-5, OFF M & P
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
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax(:,1:2:end), 'YTick', []);


%% Divergence 5, 2012-09-13-2, OFF M & Ps
piece = pieces('2012-09-13-2');
conerun = piece.d01s;
rasterrun = piece.d04cm;
map = piece.d04cm.mapd01s;
rgcs = [];
rgcs(end+1) = 7516; % offM
rgcs(end+1) = 6993; % offP
rgcs(end+1) = 7666; % offP
cones2plot = 1;
bounds = [80 150 230 290];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax, 'YLim', [0 60]);
set(rax(:,1:end-1), 'YTick', []);


%% Divergence 6, 2012-09-13-2, OFF M & P, edge
piece = pieces('2012-09-13-2');
conerun = piece.d01s;
rasterrun = piece.d04cm;
map = piece.d04cm.mapd01s;
rgcs = [];
rgcs(end+1) = 6796; % offM
rgcs(end+1) = 6798; % offP
cones2plot = [1 3 4];
bounds = [115 165 265 305];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax(:,1:2:end), 'YTick', []);


%% offM offP, but offM too small and transient
% piece = pieces('2012-09-13-2');
% conerun = piece.d01s;
% rasterrun = piece.d04cm;
% map = piece.d04cm.mapd01s;
% rgcs = [];
% rgcs(end+1) = 2117; % offM
% rgcs(end+1) = 1906; % offP
% cones2plot = 2;
% bounds = [60 115 70 120];


%% Divergence 7, 2012-09-13-2, OFF M & P
piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 2386; % offM
rgcs(end+1) = 2506; % offP
cones2plot = [2 4];
bounds = [135 190 65 120];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax, 'YLim', [0 40]);
set(rax(:,1:end-1), 'YTick', []);


%% Divergence 8, 2012-09-13-2, OFF M & Ps
piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 6991; % offM
rgcs(end+1) = 6841; % offP
rgcs(end+1) = 6601; % offP
cones2plot = 1:4;
bounds = [105 185 230 285];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax, 'YLim', [0 65]);
set(rax(:,1:end-1), 'YTick', []);


%% Divergence 9, 2012-09-13-2, OFF M & somewhat patchy P
piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 7576; % offM
rgcs(end+1) = 376; % offP
cones2plot = 1;
bounds = [45 125 195 250];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax(:,1:2:end), 'YTick', []);


%% Too patchy
% piece = pieces('2012-09-13-2');
% conerun = piece.d05s;
% rasterrun = piece.d07;
% map = piece.d07.mapd05s;
% rgcs = [];
% rgcs(end+1) = 2942; % offM
% rgcs(end+1) = 2896; % offP
% cones2plot = 1:4;
% bounds = [];
% [rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);


%% Too peripheral in the parasol
% piece = pieces('2012-09-13-2');
% conerun = piece.d05s;
% rasterrun = piece.d07;
% map = piece.d07.mapd05s;
% rgcs = [];
% rgcs(end+1) = 1396; % offM
% rgcs(end+1) = 646; % offP
% cones2plot = 1:4;
% bounds = [];
% [rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);


%% Divergence 10, 2012-09-13-2, OFF M & patchy P
piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 1966; % offM
rgcs(end+1) = 1981; % offP
cones2plot = 2:3;
bounds = [27 87 30 95];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

set(lax(:,2:end),   'YTick', []);
set(rax, 'YLim', [0 55]);
set(rax(:,1:end-1), 'YTick', []);


%% Divergence 11, 2012-09-13-2, OFF M & P, hard to get well isolated cones
piece = pieces('2012-09-13-2');
conerun = piece.d05s;
rasterrun = piece.d07;
map = piece.d07.mapd05s;
rgcs = [];
rgcs(end+1) = 2116; % offM
rgcs(end+1) = 2117; % offP
cones2plot = 3;
bounds = [63 115 72 126];
[rfax outstructs bounds] = divergence_compound_plot(conerun, rasterrun, map, rgcs, cones2plot, bounds);
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax, 'XTick', []);
set(rax(1:end-1,:), 'XTick', []);

%% Guess bounds
bounds = [Inf 0 Inf 0];
for cellid = rgcs
    cellnum = get_cell_indices(conerun, cellid); fit = conerun.stas.fits{cellnum}; mean = fit.mean; sd = fit.sd;
    xmin = mean(1) - sd(1); xmax = mean(1) + sd(1); ymin = mean(2) - sd(2); ymax = mean(2) + sd(2);
    bounds([1 3]) = min([bounds([1 3]); xmin ymin]); bounds([2 4]) = max([bounds([2 4]); xmax ymax]);
end; bounds


%% Cleanup
lax = cellfun(@(s)(s.AX(1)), outstructs);
rax = cellfun(@(s)(s.AX(2)), outstructs);
set(lax(:,2:end), 'YTick', []);
set(rax(:,1:end-1), 'YTick', []);

set(gcf, 'PaperType', 'tabloid');
orient(gcf, 'landscape');

% Match Y lims across whole fig?
% setaxesy(cellfun(@(s)(s.AX(2)), outstructs))