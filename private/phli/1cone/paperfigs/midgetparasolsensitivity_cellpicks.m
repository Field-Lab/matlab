% Basics
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
streamedopts = loadopts;
streamedopts.load_sta = true;
streamedopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
streamedopts.set_polarities = {'guess', true};
pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');
keep_vars = {'pieces'; 'loadopts'; 'streamedopts'; 'keep_vars'};

%% 2011-06-24-6
piece = '2011-06-24-6';
d.d05s = load_data([piece '/streamed/data005/data005'], streamedopts);
d.d08s = load_data([piece '/streamed/data008/data008'], streamedopts);
d.d09  = load_data([piece '/data009'], loadopts);
d.d11  = load_data([piece '/data011'], loadopts);
d.d13  = load_data([piece '/data013'], streamedopts);
d.d09.mapd05s = map_ei(d.d05s, d.d09);
d.d13.mapd05s = map_ei(d.d05s, d.d13);
d.d11.mapd08s = map_ei(d.d08s, d.d11);
d.d13.mapd08s = map_ei(d.d08s, d.d13);

% Convert old style additivity maps into format expected by cr_compound_plot
% Map 0 and 1 are the two cones, map 2 is both
d.d09 = read_stim_lisp_output(d.d09, '2011-06-24-6_data005_1and2');
    stimstruct = d.d09.stimulus;
    stimstruct.numcones = 2;
    stimstruct.mapims = {combine_maps(stimstruct.mapims(1:2))};
    stimstruct.mapnyc     = {cellfun(@(C)(C{1}), stimstruct.mapnyc(1:2),     'UniformOutput', false)};
    stimstruct.mapnycpoly = {cellfun(@(C)(C{1}), stimstruct.mapnycpoly(1:2), 'UniformOutput', false)};
    rgbs = stimstruct.rgbs;
    stimstruct.rgbs = cell(1, length(stimstruct.pulses));
    stimstruct.rgbs(stimstruct.maps == 0) = cellfun(@(A)([A; 0 0 0]), mat2cell(rgbs(stimstruct.maps == 0,:), ones(1,sum(stimstruct.maps == 0)), 3), 'UniformOutput', false);
    stimstruct.rgbs(stimstruct.maps == 1) = cellfun(@(A)([0 0 0; A]), mat2cell(rgbs(stimstruct.maps == 1,:), ones(1,sum(stimstruct.maps == 1)), 3), 'UniformOutput', false);
    stimstruct.rgbs(stimstruct.maps == 2) = cellfun(@(A)([A;     A]), mat2cell(rgbs(stimstruct.maps == 2,:), ones(1,sum(stimstruct.maps == 2)), 3), 'UniformOutput', false);
    d.d09.stimulus = stimstruct;
        d.d11 = read_stim_lisp_output(d.d11, '2011-06-24-6_f08_1and2');
        stimstruct = d.d11.stimulus;
        stimstruct.numcones = 2;
        stimstruct.mapims = {combine_maps(stimstruct.mapims(1:2))};
        stimstruct.mapnyc     = {cellfun(@(C)(C{1}), stimstruct.mapnyc(1:2),     'UniformOutput', false)};
        stimstruct.mapnycpoly = {cellfun(@(C)(C{1}), stimstruct.mapnycpoly(1:2), 'UniformOutput', false)};
        rgbs = stimstruct.rgbs;
        stimstruct.rgbs = cell(1, length(stimstruct.pulses));
        stimstruct.rgbs(stimstruct.maps == 0) = cellfun(@(A)([A; 0 0 0]), mat2cell(rgbs(stimstruct.maps == 0,:), ones(1,sum(stimstruct.maps == 0)), 3), 'UniformOutput', false);
        stimstruct.rgbs(stimstruct.maps == 1) = cellfun(@(A)([0 0 0; A]), mat2cell(rgbs(stimstruct.maps == 1,:), ones(1,sum(stimstruct.maps == 1)), 3), 'UniformOutput', false);
        stimstruct.rgbs(stimstruct.maps == 2) = cellfun(@(A)([A;     A]), mat2cell(rgbs(stimstruct.maps == 2,:), ones(1,sum(stimstruct.maps == 2)), 3), 'UniformOutput', false);
        d.d11.stimulus = stimstruct;
pieces(piece) = d;
leave(keep_vars{:});

piece = '2011-06-24-6';
sensprintpath = printpath('sens', piece);
d = pieces('2011-06-24-6');
conerun = d.d05s;
conerun.rgcs = [50 78 64 137 158 120 15 29 42]; % offM
rasterrun = d.d09;
rasterrun.rgcs = d.d09.mapd05s(get_cell_indices(conerun, conerun.rgcs));
stablerun = d.d13;
stablerun.rgcs = d.d13.mapd05s(get_cell_indices(conerun, conerun.rgcs));
[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers, 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', sensprintpath);

crprintpath = [];
d = pieces('2011-06-24-6');
conerun = d.d08s;
conerun.rgcs = [34 75 88 96 125 198]; % offP
rasterrun = d.d11;
rasterrun.rgcs = d.d11.mapd08s(get_cell_indices(conerun, conerun.rgcs));
stablerun = d.d13;
stablerun.rgcs = d.d13.mapd08s(get_cell_indices(conerun, conerun.rgcs));
[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers, 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', sensprintpath);


%% 2012-09-24-5
% These look usable, but only 2 negative contrasts...

piece = '2012-09-24-5';
d.d03s = load_data([piece '/data003'], streamedopts);
d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], streamedopts);
d.d06cm = read_stim_lisp_output(d.d06cm);
d.d03_06_07.mapd03s = map_ei(d.d03s, d.d03_06_07);
pieces(piece) = d;
leave(keep_vars{:});

piece = '2012-09-24-5';
sensprintpath = printpath('sens', piece);
d = pieces('2012-09-24-5');
conerun = d.d03s;
%conerun.rgcs = [3677 4591 2596 346 3631 3646 5206 6931 7231 3152 2851 856]; % Need decontamination
conerun.rgcs = [3677 4591 3631 3646 5206 6931 7231 856];
conerun.rgcs = [conerun.rgcs 4172 3647 887 722 2971 3616 2521];
rasterrun = d.d06cm;
rasterrun.rgcs = d.d03_06_07.mapd03s(get_cell_indices(conerun, conerun.rgcs));
stablerun = d.d07cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', sensprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d06cm d.d03s d.d07cm] = deal(cleared{:});

% Off Midgets
%   3677, cone 3
%   4591, cone 3
%   2596, cone 3; DA decontaminated; need to rerun Vision to use
%    346, cone 1; DA decontaminated; need to rerun Vision to use
%   3631, cone 4
%   3646, cone 3
%   5206, cone 4
%   6931, cone 2
%   7231, cone 1
%   3152, cone 1; DA decontaminated; need to rerun Vision to use
%   2851, cone 4; DA decontaminated; need to rerun Vision to use
%    856, cone 4
%     77, cone 4?
%   1051, cone 1?

% Off Parasols
%   4172, cone 2
%   3647, cone 4
%    887, cone 3
%    722, cone 1
%   2971, cone 1
%   3616, cone 1
%   2521, cone 4
%   4801, cone 1? 4?
%   4456, cone 2? 3? 4?
%   2941, cone 1?
%   1459, cone 1, not really good enough



%% Searching picks
d = pieces('2012-09-24-5');
conerun = d.d03s;
rasterrun = d.d06cm;
map = d.d03_06_07.mapd03s;

edge = [77 421 466 856 1051 1111 2641 2851 3152 3211 4036 4141 4336 5071 5191];
center = [7231, 6931 5206 3646, 3631 346 2596 4591 3677];

offP = [2521 2971 3616 722 4801 1459 887 2941 4456 3841 3647 4172];

%% Searching plotting

% Pick trials (same d04cm d06cm)
urgb = rasterrun.stimulus.urgb;
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
ploturgbs = [num2cell(fliplr(singleon)); num2cell(singleoff)];

% Colors to match regions
colors = jet(4);
colors = mat2cell(colors, [1 1 1 1], 3);
colors = num2cell(colors);
colors = [colors'; colors'];

% Plot
for cellid = center
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'add_width', 2, 'start', -0.2, 'color', {[0.5 0.5 0.5]}, 'hist_color', colors, 'hist_line_width', 2);

    sanesubplot(1, 3, {1 3});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end