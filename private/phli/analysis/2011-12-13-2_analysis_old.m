%% Basics
piece = '2011-12-13-2';

%% Cell picks
onM_d00  = [1561 1608 1817 1982 2283 2612 2658 2716 2822 3016 3196 3617 3961 4517 5191 6392];
offM_d00 = [3937 812 1006 707 5419];

onM_d04 = [601 1562 1681 1968 2283 2658 2821 3244 4562 6542 6601 7336];
offM_d04 = [541 5161];

offM_d08 = [541 1321 1351 2251 3031 3586 4576 5162 7486];
onM_d08 = [856];

offM_d18 = [811 916 1321 1351 2658 3586 3632 4141 4576 4861 5357 6905];
onM_d18 = [601 2661 2821 2822];


%% Data 
d00s = load_data([piece '/streamed/data000-0/data000-0']);
d04s = load_data([piece '/streamed/data004-0/data004-0']);
d08s = load_data([piece '/streamed/data008-0/data008-0']);
d14s = load_data([piece '/streamed/data014-0/data014-0']);
d18s = load_data([piece '/streamed/data018-0/data018-0']);
d23s = load_data([piece '/streamed/data023-0/data023-0']);
d03 = load_data([piece, '/data003']);
d05 = load_data([piece, '/data005']);
d07 = load_data([piece, '/data007']);
d10 = load_data([piece, '/data010']);
d12 = load_data([piece, '/data012']);
d13 = load_data([piece, '/data013']);
d20 = load_data([piece, '/data020']);
d21 = load_data([piece '/data021']);
d22 = load_data([piece '/data022']);

d00cm0 = load_data([piece '/data000_001_003_004-norefit/data000-from-data000_data001_data003_data004/data000-from-data000_data001_data003_data004']);
d01cm0 = load_data([piece '/data000_001_003_004-norefit/data001-from-data000_data001_data003_data004/data001-from-data000_data001_data003_data004']);
d03cm0 = load_data([piece '/data000_001_003_004-norefit/data003-from-data000_data001_data003_data004/data003-from-data000_data001_data003_data004']);
d04cm0 = load_data([piece '/data000_001_003_004-norefit/data004-from-data000_data001_data003_data004/data004-from-data000_data001_data003_data004']);

d04cm = load_data([piece '/data004_data007_data008-norefit/data004-from-data004_data007_data008/data004-from-data004_data007_data008']);
d07cm = load_data([piece '/data004_data007_data008-norefit/data007-from-data004_data007_data008/data007-from-data004_data007_data008']);
d08cm = load_data([piece '/data004_data007_data008-norefit/data008-from-data004_data007_data008/data008-from-data004_data007_data008']);

d08cm2 = load_data([piece '/data008_data012_data013_data014-norefit/data008-from-data008_data012_data013_data014/data008-from-data008_data012_data013_data014']);
d12cm2 = load_data([piece '/data008_data012_data013_data014-norefit/data012-from-data008_data012_data013_data014/data012-from-data008_data012_data013_data014']);
d13cm2 = load_data([piece '/data008_data012_data013_data014-norefit/data013-from-data008_data012_data013_data014/data013-from-data008_data012_data013_data014']);
d14cm2 = load_data([piece '/data008_data012_data013_data014-norefit/data014-from-data008_data012_data013_data014/data014-from-data008_data012_data013_data014']);

d18cm = load_data([piece '/data018-from-data018_data021_data022']);
d21cm = load_data([piece '/data021-from-data018_data021_data022']);
d22cm = load_data([piece '/data022-from-data018_data021_data022']);

d07da = load_data([piece, '/data007-from-data007DA']);


%% Basic load data
dataruns = {'d00s' 'd04s' 'd08s' 'd14s' 'd18s' 'd23s' 'd03' 'd05' 'd07' 'd10' 'd12' 'd13' 'd20' 'd21' 'd22' 'd00cm0' 'd04cm' 'd07cm' 'd08cm' 'd08cm2' 'd12cm2' 'd13cm2' 'd14cm2' 'd18cm' 'd21cm' 'd22cm', 'd07da'};
for i = 1:length(dataruns)
    if ~exist(dataruns{i}, 'var'), continue; end
    
    datarun = evalin('caller', dataruns{i});
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% EIs
dataruns = {'d00s' 'd04s' 'd08s' 'd14s' 'd18s' 'd03' 'd05' 'd07' 'd10' 'd12' 'd13' 'd20' 'd21' 'd22' 'd00cm0' 'd04cm' 'd08cm2' 'd18cm' 'd21cm' 'd22cm' 'd07da'};
for i = 1:length(dataruns)
    if ~exist(dataruns{i}, 'var'), continue; end

    datarun = evalin('caller', dataruns{i});
    datarun = load_ei(datarun, 'all', 'keep_java_ei', false);
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% STAs and cones
dataruns = {'d00s' 'd04s' 'd08s' 'd14s' 'd18s' 'd23s' 'd04cm' 'd08cm' 'd08cm2' 'd14cm2'};
for i = 1:length(dataruns)
    if ~exist(dataruns{i}, 'var'), continue; end

    datarun = evalin('caller', dataruns{i});
    
    datarun = load_sta(datarun,'load_sta',[], 'keep_java_sta', false, 'guess_stimulus', false);
    datarun = set_polarities(datarun);
    datarun = get_sta_fits_from_vision(datarun);
    
    assignin('caller', dataruns{i}, datarun);
end


% Normal cone loading
dataruns = {'d00s' 'd04s' 'd08s' 'd14s' 'd18s' 'd23s'};
for i = 1:length(dataruns)
    if ~exist(dataruns{i}, 'var'), continue; end

    datarun = evalin('caller', dataruns{i});

    datarun = load_cones(datarun, 1);
    
    assignin('caller', dataruns{i}, datarun);
end


% Mosaic generation
for i = 1:length(dataruns)
    if ~exist(dataruns{i}, 'var'), continue; end

    datarun = evalin('caller', dataruns{i});

    datarun = make_mosaic_struct(datarun);
    datarun = make_voronoi_masks(datarun);
    
    assignin('caller', dataruns{i}, datarun);
end

clear dataruns datarun i


%% EI map d00cm to d00s (U/D/C/R, ON midgets)
cell_list_map = map_ei(d00s, d00cm0, 'master_cell_type', {3, 4});
offM_fd04s_d00cm0 = cell_list_map(get_cell_indices(d00s, offM_d00));
onM_fd04s_d00cm0 = cell_list_map(get_cell_indices(d00s, onM_d00));
clear cell_list_map;


%% EI map d04cm to d04s (U/D/C/R, ON midgets)
cell_list_map = map_ei(d04s, d04cm, 'master_cell_type', {3, 4});
midgets_fd04s_d04cm = cell_list_map(get_cell_indices(d04s, [offM_d04 onM_d04]));
clear cell_list_map;


%% EI map d07da to d04s
cell_list_map = map_ei(d04s, d07da, 'master_cell_type', {4});
d07da.offM_fd04 = cell_list_map(get_cell_indices(d04s, offM_d04));
clear cell_list_map;

cell_list_map = map_ei(d04s, d07da, 'master_cell_type', {3}, 'electrode_threshold', 3, 'corr_threshold', 0.85);
d07da.onM_fd04  = cell_list_map(get_cell_indices(d04s, onM_d04));
clear cell_list_map;

% Refined from various map tweakings
d07da.onM_fd04 = {872 1564 [] 1516 [] [] [] 3107 [] [] 6618 []};
% notes/dupes:   {                             ?       6589   }


%% EI map d07 to d04s
cell_list_map = map_ei(d04s, d07, 'master_cell_type', {3, 4});
midgets_fd04s_d07 = cell_list_map(get_cell_indices(d04s, [offM_d04 onM_d04]));
clear cell_list_map;


%% EI map d08s to d04s
cell_list_map = map_ei(d04s, d08s, 'master_cell_type', {3, 4});
midgets_fd04s_d08s = cell_list_map(get_cell_indices(d04s, [offM_d04 onM_d04]));
clear cell_list_map;


%% EI map d08cm2 to d08s
cell_list_map = map_ei(d08s, d08cm2, 'master_cell_type', {3, 4});
offM_fd08s_d08cm2 = cell_list_map(get_cell_indices(d08s, offM_d08));
onM_fd08s_d08cm2  = cell_list_map(get_cell_indices(d08s, onM_d08));
clear cell_list_map;


%% EI map d10 to d08s (allcones)
cell_list_map = map_ei(d08s, d10, 'master_cell_type', {3, 4});
midgets_fd08s_d10 = cell_list_map(get_cell_indices(d08s, [offM_d08 onM_d08]));
clear cell_list_map;


%% EI map d12 to d08s (U/D/C/R)
cell_list_map = map_ei(d08s, d12, 'master_cell_type', {3, 4});
midgets_fd08s_d12 = cell_list_map(get_cell_indices(d08s, [offM_d08 onM_d08]));
clear cell_list_map;


%% EI map d13 to d08s (U/D/C/R)
cell_list_map = map_ei(d08s, d13, 'master_cell_type', {3, 4});
midgets_fd08s_d13 = cell_list_map(get_cell_indices(d08s, [offM_d08 onM_d08]));
clear cell_list_map;


%% EI map d14s to d08s (U/D/C/R, allcones)
cell_list_map = map_ei(d08s, d14s, 'master_cell_type', {3, 4});
midgets_fd08s_d14s = cell_list_map(get_cell_indices(d08s, [offM_d08 onM_d08]));
clear cell_list_map;

% Manual...
midgets_fd08s_d14s{2} = 1322;


%% Null stimulus
% cell_list_map = map_ei(d18s, d20, 'master_cell_type', {3, 4});
% midgets_fd18s_d20 = cell_list_map(get_cell_indices(d18s, [offM_d18 onM_d18]));
% clear cell_list_map;
% 
% % Manual...
% midgets_fd18s_d20{1} = 922;
% midgets_fd18s_d20{3} = 1321;
% midgets_fd18s_d20{5} = 2972;    % duplicates: 2971 2656
% midgets_fd18s_d20{6} = 3588;    % duplicates: 3589
% midgets_fd18s_d20{7} = 3738;    % duplicates: 3739 3741 3742
% midgets_fd18s_d20{8} = 4143;    % ???
% midgets_fd18s_d20{9} = 4576;    %
% midgets_fd18s_d20{11} = 5491;   % 
% midgets_fd18s_d20{12} = 6901;   % duplicates: 6902 6904
% 
% midgets_fd18s_d20{16} = 3036;   % ????????


%% EI map d18cm to d18s (null stimulus)
cell_list_map = map_ei(d18s, d18cm, 'master_cell_type', {3 4});
midgets_fd18s_d18cm = cell_list_map(get_cell_indices(d18s, [offM_d18 onM_d18]));
clear cell_list_map;


%% Stimuli
s03 = read_stim_lisp_output(piece, 's03');
s05 = read_stim_lisp_output(piece, 's05');
s12 = read_stim_lisp_output(piece, 's12');
s13 = read_stim_lisp_output(piece, 's13');
s20 = read_stim_lisp_output(piece, 's20');


%% Build cone flash response templates
%  Ideally these would have as many runs/cells/trials as possible, so this is just a start.
onp  = @(rgb)(all(rgb(:) >= 0) && any(rgb(:) > 0));
offp = @(rgb)(all(rgb(:) <= 0) && any(rgb(:) < 0));

runs = {d07cm d12 d13};
triggers = {d07cm.triggers(1:2:end) d12.triggers(1:2:end) d13.triggers(1:2:end)};
opts.stop = 0.8;
opts.stimstructs = {s07.rgbs s12.rgbs s13.rgbs};
off_midgets = {midgets_fd04s_d04cm(1:2) midgets_fd08s_d12(1:9) midgets_fd08s_d13(1:9)};
[templates.off_midget.on.x, templates.off_midget.on.y]  = build_response_template(runs, triggers, off_midgets, opts, 'predicates', onp);
[templates.off_midget.off.x,templates.off_midget.off.y] = build_response_template(runs, triggers, off_midgets, opts, 'predicates', offp);

runs = {d07cm d12};
triggers = {d07cm.triggers(1:2:end) d12.triggers(1:2:end)};
opts.stop = 0.8;
opts.stimstructs = {s07.rgbs s12.rgbs};
on_midgets = {midgets_fd04s_d04cm([3:10 12:14]) midgets_fd08s_d12(10)};
[templates.on_midget.on.x, templates.on_midget.on.y]  = build_response_template(runs, triggers, on_midgets, opts, 'predicates', onp);
[templates.on_midget.off.x,templates.on_midget.off.y] = build_response_template(runs, triggers, on_midgets, opts, 'predicates', offp);
clear runs on_midgets triggers opts

% Boxcar templates
templates.box.x = templates.on_midget.on.x;
templates.box.y = zeros(size(templates.box.x));
templates.box.y(templates.box.x > 0.05 & templates.box.x < 0.15) = 1;
templates.box.y = templates.box.y ./ sum(templates.box.y);


%% ww07da
udprintpath = []; %printpath('ud');
d07da = read_stim_lisp_output(d07da);

ww07da.wwrun = d07da;
ww07da.wwrun.wwrgcs = [d07da.offM_fd04 d07da.onM_fd04];
ww07da.conerun = d04s;
ww07da.conerun.wwrgcs = [offM_d04 onM_d04];
ww07da.stablerun = d08s;
ww07da.stablerun.wwrgcs = midgets_fd04s_d08s;
ww_compound_plot(ww07da, 'all', 'triggers', ww07da.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear udprintpath


%% ww07cm
udprintpath = []; %printpath('ud');
d07cm = read_stim_lisp_output(d07cm);

ww07cm.wwrun = d07cm;
ww07cm.wwrun.wwrgcs = midgets_fd04s_d04cm;
ww07cm.conerun = d04cm;
ww07cm.conerun.wwrgcs = cell2mat(midgets_fd04s_d04cm);
ww07cm.stablerun = d08cm;
ww07cm.stablerun.wwrgcs = midgets_fd04s_d04cm;
ww_compound_plot(ww07cm, 'all', 'triggers', ww07cm.wwrun.triggers(1:2:end), 'printpath', udprintpath);
%ww_compound_plot(ww07cm, [1:3 9 11:13], 'triggers', ww07cm.wwrun.triggers(1:2:end));
clear udprintpath

% Very interesting ONs:
%    3  id 601
%    9  id2821
%   11  id6541
%   12  id6601
%   13  id7336

% Okay ONs:
%    4  id1564

% Not so good ONs:
%    5  id1683
%    6  id1968
%    7  id2284
%   10  id3245


% OFFs are nice, although some stability contamination


%% Separate for figures
ww_compound_plot(ww07cm, 2, 'triggers', ww07cm.wwrun.triggers(1:2:end), 'scaled_up', 10, 'az_pad_factor', 3.5);
ww_compound_plot(ww07cm, 3, 'triggers', ww07cm.wwrun.triggers(1:2:end), 'scaled_up', 10, 'az_pad_factor', 3);


%% d07 U/D/C/R
d04s = load_data([piece '/streamed/data004-0/data004-0'], ...
                struct('load_params', true, 'load_neurons', true, 'load_ei', true, 'verbose', true));
d04cm = load_data([piece '/data004_data007_data008-norefit/data004-from-data004_data007_data008/data004-from-data004_data007_data008'], ...
                struct('load_params', true, 'load_neurons', true, 'load_ei', true, 'verbose', true));
d07cm = load_data([piece '/data004_data007_data008-norefit/data007-from-data004_data007_data008/data007-from-data004_data007_data008'], ...
                struct('load_params', true, 'load_neurons', true, 'verbose', true));
d08cm = load_data([piece '/data004_data007_data008-norefit/data008-from-data004_data007_data008/data008-from-data004_data007_data008'], ...
                struct('load_params', true, 'load_neurons', true, 'verbose', true));


% EI map d04cm to d04s (U/D/C/R, ON midgets)
cell_list_map = map_ei(d04s, d04cm, 'master_cell_type', {3});
onmidgets_fd04s_d04cm = cell_list_map(get_cell_indices(d04s, [onM_d04]));
cell_list_map = map_ei(d04s, d04cm, 'master_cell_type', {4});
offmidgets_fd04s_d04cm = cell_list_map(get_cell_indices(d04s, [offM_d04]));
clear cell_list_map;


s07 = read_stim_lisp_output(piece, 's07');
s07 = parse_stim_rgbs(s07);


% stim = 'decr';
% repr = 'incr';
% urgbs = {};
% for i = 1:s20.numcones
%     for j = 1:s20.numcones
%         if i == j
%             urgbs{i,j} = find(s20.urgb.singles & s20.urgb.(stim)(i,:));
%         else
%             urgbs{i,j} = find(s20.urgb.doubles & s20.urgb.(stim)(i,:) & s20.urgb.(repr)(j,:));
%         end
%     end
% end; clear i j stim repr singles doubles
% for i = 1:length(offmidgets_fd18s_d18cm)
%     figure;
%     urgb_raster_subplot(d20cm, offmidgets_fd18s_d18cm{i}, d20cm.triggers(1:2:end), s20, urgbs, 'hist_line_width', 1.5, 'hist_color', 'b');
% end
% clear urgbs i;


s = s07;
stim = 'incr';
repr = 'decr';
urgbs = {};
for i = 1:s.numcones
    for j = 1:s.numcones
        if i == j
            urgbs{i,j} = find(s.urgb.singles & s.urgb.(stim)(i,:));
        else
            urgbs{i,j} = find(s.urgb.doubles & s.urgb.(stim)(i,:) & s.urgb.(repr)(j,:));
        end
    end
end; clear i j stim repr s

for i = 1:length(onmidgets_fd04s_d04cm)
    cell = onmidgets_fd04s_d04cm{i};
    if isempty(cell), continue; end
    
    figure;
    urgb_raster_subplot(d07cm, cell, d07cm.triggers(1:2:end), s07, urgbs, 'hist_line_width', 1.5, 'hist_color', 'b');
end
clear i cell


%% C/R analysis 07 using both onset and offset templates
offM_d04cm_fd04s = midgets_fd04s_d04cm(1:2);
onM_d04cm_fd04s  = midgets_fd04s_d04cm(3:end);

% Calculate matched filter responses for C/R
s07.urgbs = numcellunique(s07.rgbs);
s07.onrgbs = select(s07.urgbs, onp);
s07.offrgbs = select(s07.urgbs, offp);
[d07cm.crs.on_midget.on,  d07cm.rasterhists.on_midget.on]    = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), onM_d04cm_fd04s,  s07.rgbs, templates.on_midget.on,   s07.onrgbs,  'stop', 0.8);
[d07cm.crs.on_midget.off, d07cm.rasterhists.on_midget.off]   = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), onM_d04cm_fd04s,  s07.rgbs, templates.on_midget.off,  s07.offrgbs, 'stop', 0.8);
[d07cm.crs.off_midget.on,  d07cm.rasterhists.off_midget.on]  = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), offM_d04cm_fd04s, s07.rgbs, templates.off_midget.on,  s07.onrgbs,  'stop', 0.8);
[d07cm.crs.off_midget.off, d07cm.rasterhists.off_midget.off] = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), offM_d04cm_fd04s, s07.rgbs, templates.off_midget.off, s07.offrgbs, 'stop', 0.8);

% Find the stimuli that hit single cones alone
[s07.onrgbs1conei  s07.onrgbs1cone]  = parse_single_cone_rgb_indices(s07.onrgbs);
[s07.offrgbs1conei s07.offrgbs1cone] = parse_single_cone_rgb_indices(s07.offrgbs);

% Since there are the same number of stimuli for all 4 cones in this case,
% we can repackage cleanly; for more complicated runs we would need
% additional logic to pad some cones' row to have a proper matrix.
s07.onlum1cone  = cellfun(@(a)(sum(a,2)'), s07.onrgbs1cone,  'UniformOutput', false);
s07.offlum1cone = cellfun(@(a)(sum(a,2)'), s07.offrgbs1cone, 'UniformOutput', false);
s07.onlum1cone  = cell2mat(s07.onlum1cone);
s07.offlum1cone = cell2mat(s07.offlum1cone);
s07.onlum1conei  = cell2mat(s07.onrgbs1conei);
s07.offlum1conei = cell2mat(s07.offrgbs1conei);
s07.crsx = [s07.offlum1cone s07.onlum1cone];

for i = 1:length(onM_d04cm_fd04s)
    rgc = onM_d04cm_fd04s{i};
    if isempty(rgc), continue; end
    
    offcrs = zeros(size(s07.offlum1cone));
    for j = 1:numel(s07.offlum1conei)
        offcrs(j) = d07cm.crs.on_midget.off(i,s07.offlum1conei(j));
    end
    
    oncrs = zeros(size(s07.onlum1cone));
    for j = 1:numel(s07.onlum1conei)
        oncrs(j) = d07cm.crs.on_midget.on(i,s07.onlum1conei(j));
    end
    
    crs = [offcrs oncrs];
    [p resnorm] = normcdfxscale(crs, s07.crsx, 'plot', true);
    xscalesOnM07{i} = p(1:4) ./ p(1);
    xscaleErrOnM07(i) = resnorm;
    waitfor(gca);
end

for i = 1:length(offM_d04cm_fd04s)
    rgc = offM_d04cm_fd04s{i};
    if isempty(rgc), continue; end

    offcrs = zeros(size(s07.offlum1cone));
    for j = 1:numel(s07.offlum1conei)
        offcrs(j) = d07cm.crs.off_midget.off(i,s07.offlum1conei(j));
    end
    
    oncrs = zeros(size(s07.onlum1cone));
    for j = 1:numel(s07.onlum1conei)
        oncrs(j) = d07cm.crs.off_midget.on(i,s07.onlum1conei(j));
    end
    
    crs = [offcrs oncrs];
    [p resnorm] = normcdfxscale(crs, s07.crsx, 'plot', true);
    xscalesOffM07{i} = p(1:4) ./ p(1);
    xscaleErrOffM07(i) = resnorm;
    waitfor(gca);
end;
clear i j rgc offcrs oncrs crs p resnorm


%% C/R analysis 07 box template only
offM_d04cm_fd04s = midgets_fd04s_d04cm(1:2);
onM_d04cm_fd04s  = midgets_fd04s_d04cm(3:end);

% Calculate matched filter responses for C/R
s07.urgbs = numcellunique(s07.rgbs);
s07.onrgbs = select(s07.urgbs, onp);
s07.offrgbs = select(s07.urgbs, offp);
s07.blankrgb = select(s07.urgbs, @(a) (all(a(:) == 0)));
[d07cm.crs.on_midget.on,   d07cm.rasterhists.on_midget.on]   = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), onM_d04cm_fd04s,  s07.rgbs, templates.box, s07.onrgbs,  'stop', 0.8);
[d07cm.crs.on_midget.off,  d07cm.rasterhists.on_midget.off]  = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), onM_d04cm_fd04s,  s07.rgbs, templates.box, s07.offrgbs, 'stop', 0.8);
[d07cm.crs.off_midget.on,  d07cm.rasterhists.off_midget.on]  = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), offM_d04cm_fd04s, s07.rgbs, templates.box, s07.onrgbs,  'stop', 0.8);
[d07cm.crs.off_midget.off, d07cm.rasterhists.off_midget.off] = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), offM_d04cm_fd04s, s07.rgbs, templates.box, s07.offrgbs, 'stop', 0.8);
[d07cm.crs.on_midget.blank,  d07cm.rasterhists.on_midget.blank]  = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), onM_d04cm_fd04s,  s07.rgbs, templates.box,   s07.blankrgb,  'stop', 0.8);
[d07cm.crs.off_midget.blank, d07cm.rasterhists.off_midget.blank] = single_map_contrast_response(d07cm, d07cm.triggers(1:2:end), offM_d04cm_fd04s, s07.rgbs, templates.box, s07.blankrgb,  'stop', 0.8);

% Find the stimuli that hit single cones alone
[s07.onrgbs1conei  s07.onrgbs1cone]  = parse_single_cone_rgb_indices(s07.onrgbs);
[s07.offrgbs1conei s07.offrgbs1cone] = parse_single_cone_rgb_indices(s07.offrgbs);

% Since there are the same number of stimuli for all 4 cones in this case,
% we can repackage cleanly; for more complicated runs we would need
% additional logic to pad some cones' row to have a proper matrix.
s07.onlum1cone  = cellfun(@(a)(sum(a,2)'), s07.onrgbs1cone,  'UniformOutput', false);
s07.offlum1cone = cellfun(@(a)(sum(a,2)'), s07.offrgbs1cone, 'UniformOutput', false);
s07.onlum1cone  = cell2mat(s07.onlum1cone);
s07.offlum1cone = cell2mat(s07.offlum1cone);
s07.onlum1conei  = cell2mat(s07.onrgbs1conei);
s07.offlum1conei = cell2mat(s07.offrgbs1conei);
s07.crsx = [s07.offlum1cone s07.onlum1cone];

for i = 1:length(onM_d04cm_fd04s)
    rgc = onM_d04cm_fd04s{i};
    if isempty(rgc), continue; end
    
    offcrs = zeros(size(s07.offlum1cone));
    for j = 1:numel(s07.offlum1conei)
        offcrs(j) = d07cm.crs.on_midget.off(i,s07.offlum1conei(j));
    end
    
    oncrs = zeros(size(s07.onlum1cone));
    for j = 1:numel(s07.onlum1conei)
        oncrs(j) = d07cm.crs.on_midget.on(i,s07.onlum1conei(j));
    end
    
    crs = [offcrs oncrs];
    [p resnorm] = normcdfxscale(crs, s07.crsx, 'plot', true, 'colors', [1 0 0; 0 0.5 0; 0 0 1; 1 0 1]);
    xscalesOnM07{i} = p(1:4) ./ p(1);
    xscaleErrOnM07(i) = resnorm;
    hold on;
    plot(0, d07cm.crs.on_midget.blank(i), 'k.', 'MarkerSize', 15);
    waitfor(gca);
end

for i = 1:length(offM_d04cm_fd04s)
    rgc = offM_d04cm_fd04s{i};
    if isempty(rgc), continue; end

    offcrs = zeros(size(s07.offlum1cone));
    for j = 1:numel(s07.offlum1conei)
        offcrs(j) = d07cm.crs.off_midget.off(i,s07.offlum1conei(j));
    end
    
    oncrs = zeros(size(s07.onlum1cone));
    for j = 1:numel(s07.onlum1conei)
        oncrs(j) = d07cm.crs.off_midget.on(i,s07.onlum1conei(j));
    end
    
    crs = [offcrs oncrs];
    [p resnorm] = normcdfxscale(crs, s07.crsx, 'plot', true, 'colors', [1 0 0; 0 0.5 0; 0 0 1; 1 0 1]);
    xscalesOffM07{i} = p(1:4) ./ p(1);
    xscaleErrOffM07(i) = resnorm;
    hold on;
    plot(0, d07cm.crs.off_midget.blank(i), 'k.', 'MarkerSize', 15);
    waitfor(gca);
end;
clear i j rgc offcrs oncrs crs p resnorm


%%
% ww07.wwrun = d07;
% ww07.wwrun.wwrgcs = midgets_fd04s_d07;
% ww07.wwstim = s07;
% ww07.conerun = d04s;
% ww07.conerun.wwrgcs = [offM_d04 onM_d04];
% ww07.stablerun = d08s;
% ww07.stablerun.wwrgcs = midgets_fd04s_d08s;
% ww_compound_plot(ww07, 3, 'triggers', d07.triggers(1:2:end));

% Interesting ones!  See if we can get more cells mapped manually
% See concatenated mapping above


%% data010 Allcones analysis
allconesprintpath = []; % printpath('allcones');

conerun = d08s;
d10 = read_stim_lisp_output(d10, ':2011-12-13-2_f08_allcones');
d10.stimulus = parse_stim_rgbs(d10.stimulus);
d10.stimulus.triggers = d10.triggers(1:2:end);
datarun = d10;
stablerun = d14s;

conerun.rgcs = offM_d08;
datarun.rgcs = midgets_fd08s_d10;
stablerun.rgcs = midgets_fd08s_d14s;

%allcones_plot(datarun, conerun, 'printpath', allconesprintpath, 'urgbs', d10.stimulus.urgbs(4), 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
allcones_plot(datarun, conerun, 'urgbs', d10.stimulus.urgbs(4), 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
clear conerun datarun allconesprintpath


%%
udprintpath = []; %printpath('ud');
ww12cm2.wwrun = d12cm2;
ww12cm2.wwrun.wwrgcs = [offM_fd08s_d08cm2 onM_fd08s_d08cm2];
ww12cm2.wwrun.stimulus = s12;
ww12cm2.conerun = d08cm2;
ww12cm2.conerun.wwrgcs = cell2mat([offM_fd08s_d08cm2 onM_fd08s_d08cm2]);
ww12cm2.stablerun = d14cm2;
ww12cm2.stablerun.wwrgcs = [offM_fd08s_d08cm2 onM_fd08s_d08cm2];
ww_compound_plot(ww12cm2, 'all', 'triggers', d12cm2.triggers(1:2:end), 'printpath', udprintpath);
clear udprintpath


%% Separate for figures
ww_compound_plot(ww12cm2, 5, 'triggers', d12cm2.triggers(1:2:end), 'set_ylims', [0 40], 'MarkerSize', 15, 'scaled_up', 10, 'az_pad_factor', 3.5);
ww_compound_plot(ww12cm2, 9, 'triggers', d12cm2.triggers(1:2:end), 'set_ylims', [], 'MarkerSize', 15, 'scaled_up', 10, 'az_pad_factor', 3.5);



%%
% ww12.wwrun = d12;
% ww12.wwrun.wwrgcs = midgets_fd08s_d12;
% ww12.wwstim = s12;
% ww12.conerun = d08s;
% ww12.conerun.wwrgcs = [offM_d08 onM_d08];
% ww12.stablerun = d14s;
% ww12.stablerun.wwrgcs = midgets_fd08s_d14s;
% ww_compound_plot(ww12, 'all', 'triggers', d12.triggers(1:2:end));

% 10 ON is interesting, id 856 which is d04s id 858

%%
s12.urgbs = numcellunique(s12.rgbs);
s12.onrgbs = select(s12.urgbs, onp);
s12.offrgbs = select(s12.urgbs, offp);

%% Calculate matched filter responses for C/R
s12.urgbs = numcellunique(s12.rgbs);
s12.onrgbs = select(s12.urgbs, onp);
s12.offrgbs = select(s12.urgbs, offp);
[d12.crs.off_midget.on,  d12.rasterhists.off_midget.on]  = single_map_contrast_response(d12, d12.triggers(1:2:end), midgets_fd08s_d12(1:9), s12.rgbs, templates.off_midget.on,  s12.onrgbs,  'stop', 0.8);
[d12.crs.off_midget.off, d12.rasterhists.off_midget.off] = single_map_contrast_response(d12, d12.triggers(1:2:end), midgets_fd08s_d12(1:9), s12.rgbs, templates.off_midget.off, s12.offrgbs, 'stop', 0.8);

s13.urgbs = numcellunique(s13.rgbs);
s13.offrgbs = select(s13.urgbs, offp);
[d13.crs.off_midget.off, d13.rasterhists.off_midget.off] = single_map_contrast_response(d13, d13.triggers(1:2:end), midgets_fd08s_d13(1:9), s13.rgbs, templates.off_midget.off, s13.offrgbs, 'stop', 0.8);


%% d20 ww
udprintpath = []; %printpath('ud');

% EI map d20 to d18s
cell_list_map = map_ei(d18s, d20, 'master_cell_type', {3 4});
midgets_fd18s_d20 = cell_list_map(get_cell_indices(d18s, [offM_d18 onM_d18]));
clear cell_list_map;

% EI map d23s to d18s
% cell_list_map = map_ei(d18s, d23s, 'master_cell_type', {3 4});
% midgets_fd18s_d23s = cell_list_map(get_cell_indices(d18s, [offM_d18 onM_d18]));
% clear cell_list_map;
d20 = read_stim_lisp_output(d20);


ww20.wwrun = d20;
ww20.wwrun.wwrgcs = midgets_fd18s_d20;
ww20.conerun = d18s;
ww20.conerun.wwrgcs = [offM_d18 onM_d18];
ww20.stablerun = d23s;
ww_compound_plot(ww20, 'all', 'triggers', d20.triggers(1:2:end), 'printpath', udprintpath);
clear udprintpath

% ?? Not sure where this note came from...
% 5357 pretty interesting...


%% d20cm ww
udprintpath = []; %printpath('ud');

d18s = load_data([piece '/streamed/data018-0/data018-0'], ...
                struct('load_sta', true, 'load_params', true, 'load_neurons', true, 'load_ei', true, 'verbose', true));
d18cm = load_data([piece '/data001-018-020-023-norefit/data018-from-data001_data018_data020_data023/data018-from-data001_data018_data020_data023'], ...
                struct('load_params', true, 'load_neurons', true, 'load_ei', true, 'verbose', true));
d20cm = load_data([piece '/data001-018-020-023-norefit/data020-from-data001_data018_data020_data023/data020-from-data001_data018_data020_data023'], ...
                struct('load_params', true, 'load_neurons', true, 'verbose', true));
d23cm = load_data([piece '/data001-018-020-023-norefit/data023-from-data001_data018_data020_data023/data023-from-data001_data018_data020_data023'], ...
                struct('load_params', true, 'load_neurons', true, 'verbose', true));

% EI map d18cm to d18s
cell_list_map = map_ei(d18s, d18cm, 'master_cell_type', {4});
offmidgets_fd18s_d18cm = cell_list_map(get_cell_indices(d18s, offM_d18));
clear cell_list_map;
offmidgets_fd18s_d18cm{3} = 1321;
offmidgets_fd18s_d18cm{6} = 3586;
offmidgets_fd18s_d18cm{7} = 3737;
offmidgets_fd18s_d18cm{10} = 4862;
offmidgets_fd18s_d18cm{11} = 5491;
offmidgets_fd18s_d18cm{12} = 6901;

d18s = set_polarities(d18s);
d18s = get_sta_fits_from_vision(d18s);
d20cm = read_stim_lisp_output(d20cm);

d23cm = set_polarities(d23cm);
d23cm = load_sta(d23cm, 'load_sta', []);
d23cm = get_sta_fits_from_vision(d23cm);

ww20cm.wwrun = d20cm;
ww20cm.wwrun.wwrgcs = offmidgets_fd18s_d18cm;
ww20cm.conerun = d18s;
ww20cm.conerun.wwrgcs = [offM_d18];
ww20cm.stablerun = d23cm;
ww20cm.stablerun.wwrgcs = offmidgets_fd18s_d18cm;
ww_compound_plot(ww20cm, 'all', 'triggers', d20cm.triggers(1:2:end), 'printpath', udprintpath);
clear udprintpath


%% d20 U/D/C/R
d18s = load_data([piece '/streamed/data018-0/data018-0'], ...
                struct('load_params', true, 'load_neurons', true, 'load_ei', true, 'verbose', true));
d18cm = load_data([piece '/data001-018-020-023-norefit/data018-from-data001_data018_data020_data023/data018-from-data001_data018_data020_data023'], ...
                struct('load_params', true, 'load_neurons', true, 'load_ei', true, 'verbose', true));
d20cm = load_data([piece '/data001-018-020-023-norefit/data020-from-data001_data018_data020_data023/data020-from-data001_data018_data020_data023'], ...
                struct('load_params', true, 'load_neurons', true, 'verbose', true));

% EI map d18cm to d18s
cell_list_map = map_ei(d18s, d18cm, 'master_cell_type', {4});
offmidgets_fd18s_d18cm = cell_list_map(get_cell_indices(d18s, offM_d18));
clear cell_list_map;
offmidgets_fd18s_d18cm{3} = 1321;
offmidgets_fd18s_d18cm{6} = 3586;
offmidgets_fd18s_d18cm{7} = 3737;
offmidgets_fd18s_d18cm{10} = 4862;
offmidgets_fd18s_d18cm{11} = 5491;
offmidgets_fd18s_d18cm{12} = 6901;

cell_list_map = map_ei(d18s, d18cm, 'master_cell_type', {3});
onmidgets_fd18s_d18cm = cell_list_map(get_cell_indices(d18s, onM_d18));
clear cell_list_map;


s20 = read_stim_lisp_output(piece, 's20');
s20 = parse_stim_rgbs(s20);


stim = 'decr';
repr = 'incr';
urgbs = {};
for i = 1:s20.numcones
    for j = 1:s20.numcones
        if i == j
            urgbs{i,j} = find(s20.urgb.singles & s20.urgb.(stim)(i,:));
        else
            urgbs{i,j} = find(s20.urgb.doubles & s20.urgb.(stim)(i,:) & s20.urgb.(repr)(j,:));
        end
    end
end; clear i j stim repr singles doubles
for i = 1:length(offmidgets_fd18s_d18cm)
    figure;
    urgb_raster_subplot(d20cm, offmidgets_fd18s_d18cm{i}, d20cm.triggers(1:2:end), s20, urgbs, 'hist_line_width', 1.5, 'hist_color', 'b');
end
clear urgbs i;


stim = 'incr';
repr = 'decr';
urgbs = {};
for i = 1:s20.numcones
    for j = 1:s20.numcones
        if i == j
            urgbs{i,j} = find(s20.urgb.singles & s20.urgb.(stim)(i,:));
        else
            urgbs{i,j} = find(s20.urgb.doubles & s20.urgb.(stim)(i,:) & s20.urgb.(repr)(j,:));
        end
    end
end; clear i j stim repr singles doubles
for i = 1:length(onmidgets_fd18s_d18cm)
    figure;
    urgb_raster_subplot(d20cm, onmidgets_fd18s_d18cm{i}, d20cm.triggers(1:2:end), s20, urgbs, 'hist_line_width', 1.5, 'hist_color', 'b');
end
clear i urgbs


%% Null stimulus analysis
for i = 1:length(offM_d18)
    if isempty(midgets_fd18s_d18cm{i}), continue; end
    
    figure;
    sanesubplot(2, 1, [1 1]);
    [res_on,    ax_on]  = rasterphli(d21cm, midgets_fd18s_d18cm{i}, d21cm.triggers(1:2:end), 'hist', true, 'hist_line_width', 2, 'MarkerSize', 5, 'plot', true);
    title(['Off Midget ' num2str(midgets_fd18s_d18cm{i})]);
    
    sanesubplot(2, 1, [2 1]);
    [res_null, ax_null] = rasterphli(d22cm, midgets_fd18s_d18cm{i}, d22cm.triggers(1:2:end), 'hist', true, 'hist_line_width', 2, 'MarkerSize', 5, 'plot', true);
    
    linkaxes([ax_on(2) ax_null(2)]);
    set(ax_on(2), 'YTickMode', 'auto');
    set(ax_null(2), 'YTickMode', 'auto');
end; clear i

clear res_on res_null ax_on ax_null


%% Redo null stimulus for figure
rgcid = midgets_fd18s_d18cm{1};
figure;

sanesubplot(2, 2, [1 2]);
[res_on,    ax_on]  = rasterphli(d21cm, rgcid, d21cm.triggers(1:2:end), 'hist', true, 'hist_line_width', 2, 'MarkerSize', 10, 'plot', true);
title(['Off Midget ' num2str(rgcid)]);

sanesubplot(2, 2, [2 2]);
[res_null, ax_null] = rasterphli(d22cm, rgcid, d22cm.triggers(1:2:end), 'hist', true, 'hist_line_width', 2, 'MarkerSize', 10, 'plot', true);
linkaxes([ax_on(2) ax_null(2)]);
set(ax_on(2), 'YTickMode', 'auto');
set(ax_null(2), 'YTickMode', 'auto');

map21 = dlmread([server_data_path piece '/2011-12-13-2_f18_nulltest/map-0001.txt']);
map22 = dlmread([server_data_path piece '/2011-12-13-2_f18_nulltest/map-0000.txt']);
nyc21 = map2manhattan(map21);
nyc22 = map2manhattan(map22);

d18cm = load_sta(d18cm, 'load_sta', rgcid);
d18cm = set_polarities(d18cm);
d18cm = get_sta_fits_from_vision(d18cm);

sanesubplot(2, 2, [1 1]);
plot_rf(d18cm, rgcid, 'fit', false, 'scale', 10, 'title', false);
hold on;
% fillmap(map21, 'colors', 'r', 'scale', [1 1]./3, 'EdgeColor', 'r');
patchmanhattan(nyc21, 'colors', [0 0 0], 'scale', [1 1]./3);
autozoom_to_fit(d18cm, rgcid, 5.4, 1, 1.2);
set(gca, 'XTick', [], 'YTick', []);

sanesubplot(2, 2, [2 1]);
plot_rf(d18cm, rgcid, 'fit', false, 'scale', 10, 'title', false);
hold on;
% fillmap(map22, 'colors', 'r', 'scale', [1 1]./3, 'EdgeColor', 'r');
patchmanhattan(nyc22, 'colors', [0 0 0], 'scale', [1 1]./3);
autozoom_to_fit(d18cm, rgcid, 5.4, 1, 1.2);
set(gca, 'XTick', [], 'YTick', []);

% Including FaceAlpha causes Matlab to export like a retard, so have to set
% transparency in postprocessing

% EdgeColor none gives you disjoint squares in Intaglio.  BARF!!!!!


%% Cone projective field; searching

ww_raster_plot(d07cm, s07, num2cell(d04cm.cell_types{3}.cell_ids), [d04cm.stas.polarities{get_cell_indices(d04cm, {3})}], 'triggers', d07cm.triggers(1:2:end), 'figure_by_cell_id', true);
% Polarities a bit strange
on_midgets_incidental_d07cm = [529, 707, 6406, 6436];

cpf07cm.wwrun = d07cm;
cpf07cm.wwrun.wwrgcs = num2cell(on_midgets_incidental_d07cm);
cpf07cm.wwstim = s07;
cpf07cm.conerun = d04cm;
cpf07cm.conerun.wwrgcs = on_midgets_incidental_d07cm;
cpf07cm.stablerun = d08cm;
cpf07cm.stablerun.wwrgcs = num2cell(on_midgets_incidental_d07cm);
ww_compound_plot(cpf07cm, 'all', 'triggers', d07cm.triggers(1:2:end));
% 707 okay, but movement
% 6406, 6436 interesting; surround?


ww_raster_plot(d07cm, s07, num2cell(d04cm.cell_types{4}.cell_ids), [d04cm.stas.polarities{get_cell_indices(d04cm, {4})}], 'triggers', d07cm.triggers(1:2:end), 'figure_by_cell_id', true);
off_midgets_incidental_d07cm = [588, 589, 811, 1351, 1591, 1742, 1906, 2191, 2656, 2839, 2971, 3243, 3406, 3751, 4621, 6260, 6455];

cpf07cm.wwrun = d07cm;
cpf07cm.wwrun.wwrgcs = num2cell(off_midgets_incidental_d07cm);
cpf07cm.wwstim = s07;
cpf07cm.conerun = d04cm;
cpf07cm.conerun.wwrgcs = off_midgets_incidental_d07cm;
cpf07cm.stablerun = d08cm;
cpf07cm.stablerun.wwrgcs = num2cell(off_midgets_incidental_d07cm);
ww_compound_plot(cpf07cm, 'all', 'triggers', d07cm.triggers(1:2:end));
nice_off_midgets_incidental_d07cm = [1351 1591 2191 2656 2971 3406 3751 4621];

ww_raster_plot(d07cm, s07, num2cell(d04cm.cell_types{2}.cell_ids), [d04cm.stas.polarities{get_cell_indices(d04cm, {2})}], 'triggers', d07cm.triggers(1:2:end), 'figure_by_cell_id', true);
off_parasols_incidental_d07cm = [466 1126 2101 2312 2926 3438 4291 4966 5071 5446 6391 6603 6752 7201];

cpf07cm.wwrun = d07cm;
cpf07cm.wwrun.wwrgcs = num2cell(off_parasols_incidental_d07cm);
cpf07cm.wwstim = s07;
cpf07cm.conerun = d04cm;
cpf07cm.conerun.wwrgcs = off_parasols_incidental_d07cm;
cpf07cm.stablerun = d08cm;
cpf07cm.stablerun.wwrgcs = num2cell(off_parasols_incidental_d07cm);
ww_compound_plot(cpf07cm, 'all', 'triggers', d07cm.triggers(1:2:end));
nice_off_parasols_incidental_d07cm = [4291];
okay_off_parasols_incidental_d07cm = [4966, 5071, 6752];
tricky_off_parasols_incidental_d07cm = [2312, 3438, 6391, 7201];
doubled_off_parasols_incidental_d07cm = [466, 6603];


ww_raster_plot(d07cm, s07, num2cell(d04cm.cell_types{9}.cell_ids), 1, 'triggers', d07cm.triggers(1:2:end), 'figure_by_cell_id', true);
% 1728, 2222, 2267, 2507, 3061, 3409, 3801, 4741, 5029, 5978, 6286, 6454, 6482, 6484, 6587, 6737, 7116
% maybe one after 7247?

% plot_rf_summaries(d04s, [offM_d04 onM_d04], 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b')
% plot_rf_summaries(d04cm, on_midgets_incidental_d07cm, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
% plot_rf_summaries(d04cm, off_midgets_incidental_d07cm, 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'g')


%% Stability summary
piece = '2011-12-13-2';
d0 = load_data(sprintf('%s/streamed/data%.3d-0/data%.3d-0', piece, 0, 0));
d0 = load_globals(d0);
d0time = d0.globals.readRDH512().getTime();

coneruns = [0 4 8 11 14 18 23];
numplots = length(coneruns) - 1;

fig = figure();
for i = 1:numplots
    subplot(2, ceil(numplots/2), i);
    
    d1run = coneruns(i);
    d2run = coneruns(i+1);
    
    d1 = load_data(sprintf('%s/streamed/data%.3d-0/data%.3d-0', piece, d1run, d1run));
    d1 = load_params(d1);
    d1 = load_cones(d1, 'cquisition');
    d1 = load_globals(d1);
    d1 = stimulus_from_globals(d1);
    
    d2 = load_data(sprintf('%s/streamed/data%.3d-0/data%.3d-0', piece, d2run, d2run));
    d2 = load_params(d2);
    d2 = load_cones(d2, 'cquisition');
    d2 = load_globals(d2);
    d2 = stimulus_from_globals(d2);
    
    overlay_cone_mosaics(d1, d2);
    axis equal tight
   
    secdiff = d2.globals.readRDH512.getTime - d1.globals.readRDH512.getTime;
    startdiff = d2.globals.readRDH512.getTime - d0time;
    title(sprintf('data%.3d vs data%.3d\n%s apart\n%s since data000', d1run, d2run, printtime(secdiff), printtime(startdiff)));
end

% Print
set(fig, 'PaperType', 'B');
orient(fig, 'landscape');
set(fig, 'Color', 'w');
set(fig,'InvertHardcopy','off')
print('-dpdf', '-painters', '-noui', printpath('stability', sprintf('stability-%s', piece)));