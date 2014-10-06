%% Basics
piece = '2011-12-13-2';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%% Streamed runs and cell picks

d.d00s = load_data([piece '/streamed/data000-0/data000-0'], staopts);
d.d00s.onM  = [1561 1608 1817 1982 2283 2612 2658 2716 2822 3016 3196 3617 3961 4517 5191 6392];
d.d00s.offM = [3937 812 1006 707 5419];

d.d04s = load_data([piece '/streamed/data004-0/data004-0'], staopts);
d.d04s.onM  = [601 1562 1681 1968 2283 2658 2821 3244 4562 6542 6601 7336];
d.d04s.offM = [541 5161];
d.d04s = load_cones(d.d04s, 'bayes');

d.d08s = load_data([piece '/streamed/data008-0/data008-0'], staopts);
d.d08s.offM = [541 1321 1351 2251 3031 3586 4576 5162 7486];
d.d08s.onM  = [856];
d.d08s = load_cones(d.d08s, 'bayes');

d.d18s = load_data([piece '/streamed/data018-0/data018-0'], staopts);
d.d18s.offM = [811 916 1321 1351 2658 3586 3632 4141 4576 4861 5357 6905];
d.d18s.onM  = [601 2661 2821 2822];


%% Good RF ON Midgets
d.d04s.rgcs = [601 2778 3137 6542 7336];
d.d00s.mapd04s = map_ei(d.d04s, d.d00s);
d.d08s.mapd04s = map_ei(d.d04s, d.d08s);

d.d11s = load_data([piece '/streamed/data011-0/data011-0'], staopts);
d.d11s.mapd04s = map_ei(d.d04s, d.d11s);

mappedruns = {d.d00s d.d08s d.d11s};
mappedrunindices = {[] [] [1] [2 3] [2 3]};
for i = 1:length(d.d04s.rgcs)
    id04 = d.d04s.rgcs(i);
    cellnum = get_cell_indices(d.d04s, id04);
    compoundrf = get_rf(d.d04s, id04);
    
    for mappedrunindex = mappedrunindices{i}
        mappedrun = mappedruns{mappedrunindex};
        mappedid = mappedrun.mapd04s{cellnum};
        if isempty(mappedid), continue; end
        compoundrf = compoundrf + get_rf(mappedrun, mappedid);
    end
    
    figure;
    dummyrun = d.d04s;
    dummyrun.stas.rfs{cellnum} = compoundrf;
    plot_rf(dummyrun, id04, 'fit', false, 'autozoom', true, 'scale', 10);
end

% Used for searching...
% mappedruns = {d.d00s d.d08s d.d11s};
% for id04 = d.d04s.rgcs(5)
%     cellnum = get_cell_indices(d.d04s, id04);
%     compoundrf = get_rf(d.d04s, id04);
% 
%     figure;
%     subplot(1, length(mappedruns)+2 , 1);
%     plot_voronoi_over_rf(d.d04s, id04, 'az_aspect_ratio', 1);
% 
%     
%     for i = 1:length(mappedruns)
%         mappedrun = mappedruns{i};
%         mappedid = mappedrun.mapd04s(cellnum);
%         mappedid = mappedid{1};
%         if isempty(mappedid), continue; end
%         
%         compoundrf = compoundrf + get_rf(mappedrun, mappedid);
%         subplot(1, length(mappedruns)+2, i+1);
%         plot_voronoi_over_rf(mappedrun, mappedid, 'voronoirun', d.d04s, 'az_aspect_ratio', 1), 
%     end
%     
%     subplot(1,length(mappedruns)+2, length(mappedruns)+2);
%     dummyrun = d.d04s;
%     dummyrun.stas.rfs{cellnum} = compoundrf;
%     plot_rf(dummyrun, id04, 'fit', false, 'autozoom', true, 'scale', 10);
% end


%% Good RF OFF Midgets
d.d04s.rgcs = [2971 3136 3586 3736];

%% U/D data007cm
udprintpath = [];
d.d07cm = load_data([piece '/data004_data007_data008-norefit/data007-from-data004_data007_data008/data007-from-data004_data007_data008'], loadopts);
d.d08cm = load_data([piece '/data004_data007_data008-norefit/data008-from-data004_data007_data008/data008-from-data004_data007_data008'], staopts);

% Map from d04s to concatenated (valid also for individuals from concatenated since norefit)
d.d04_07_08 = load_data([piece '/data004_data007_data008-norefit/data004_data007_data008-norefit'], loadopts);
cell_list_map = map_ei(d.d04s, d.d04_07_08);
d.d07cm.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offM));
d.d07cm.onM_fd04s  = cell_list_map(get_cell_indices(d.d04s, d.d04s.onM));
d.d08cm.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offM));
d.d08cm.onM_fd04s  = cell_list_map(get_cell_indices(d.d04s, d.d04s.onM));
clear cell_list_map;

% Read stimulus
d.d07cm = read_stim_lisp_output(d.d07cm);
d.d07cm.stimulus = parse_stim_rgbs(d.d07cm.stimulus);

conerun = d.d04s;
conerun.wwrgcs = [conerun.onM conerun.offM];

stablerun = d.d08cm;
stablerun.wwrgcs = [stablerun.onM_fd04s stablerun.offM_fd04s];

ud07cm.wwrun = d.d07cm;
ud07cm.wwrun.wwrgcs = [ud07cm.wwrun.onM_fd04s ud07cm.wwrun.offM_fd04s];
ud07cm.conerun = conerun;
ud07cm.stablerun = stablerun;
ww_compound_plot(ud07cm, 'all', 'triggers', ud07cm.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear conerun stablerun udprintpath

% For J Freeman subunits analysis
save_udraster_data(ud07cm, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


%% C/R data007cm
crprintpath = [];%printpath('cr', piece);

if ~isfield(d, 'd07cm'), d.d07cm = load_data([piece '/data004_data007_data008-norefit/data007-from-data004_data007_data008/data007-from-data004_data007_data008'], loadopts); end
if ~isfield(d, 'd08cm'), d.d08cm = load_data([piece '/data004_data007_data008-norefit/data008-from-data004_data007_data008/data008-from-data004_data007_data008'], loadopts); end

% Map from d04s to concatenated (valid also for individuals from concatenated since norefit)
if ~isfield(d.d07cm, 'offM_fd04s')
    if ~isfield(d, 'd04_07_08'), d.d04_07_08 = load_data([piece '/data004_data007_data008-norefit/data004_data007_data008-norefit'], loadopts); end

    cell_list_map = map_ei(d.d04s, d.d04_07_08);
    d.d07cm.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offM));
    d.d07cm.onM_fd04s  = cell_list_map(get_cell_indices(d.d04s, d.d04s.onM));
    d.d08cm.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offM));
    d.d08cm.onM_fd04s  = cell_list_map(get_cell_indices(d.d04s, d.d04s.onM));
    clear cell_list_map;
end

rasterrun = d.d07cm;
rasterrun.rgcs = [d.d07cm.offM_fd04s d.d07cm.onM_fd04s];

conerun = d.d04s;
conerun.rgcs = [conerun.offM conerun.onM];

stablerun = d.d08cm;
stablerun.rgcs = [stablerun.offM_fd04s stablerun.onM_fd04s];

% The rest of the ONs are sketchy or trivial looking...
[rasterrun conerun stablerun crs crsx ps resnorms residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d07cm d.d04s d.d08cm] = deal(cleared{:});


%% C/R data007
crprintpath = [];%printpath('cr', piece);

if ~isfield(d, 'd07'),   d.d07 = load_data([piece '/data007'], loadopts); end
if ~isfield(d, 'd08cm'), d.d08cm = load_data([piece '/data004_data007_data008-norefit/data008-from-data004_data007_data008/data008-from-data004_data007_data008'], loadopts); end

% Map from d04s to d07
if ~isfield(d.d07, 'offM_fd04s')
    cell_list_map = map_ei(d.d04s, d.d07);
    d.d07.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offM));
    d.d07.onM_fd04s  = cell_list_map(get_cell_indices(d.d04s, d.d04s.onM));
    clear cell_list_map;
end

% Map from d04s to concatenated (valid also for individuals from concatenated since norefit)
d.d04_07_08 = load_data([piece '/data004_data007_data008-norefit/data004_data007_data008-norefit'], loadopts);
cell_list_map = map_ei(d.d04s, d.d04_07_08);
d.d08cm.offM_fd04s = cell_list_map(get_cell_indices(d.d04s, d.d04s.offM));
d.d08cm.onM_fd04s  = cell_list_map(get_cell_indices(d.d04s, d.d04s.onM));
clear cell_list_map;

rasterrun = d.d07;
rasterrun.rgcs = [rasterrun.offM_fd04s rasterrun.onM_fd04s];

conerun = d.d04s;
conerun.rgcs = [conerun.offM conerun.onM];

stablerun = d.d08cm;
stablerun.rgcs = [stablerun.offM_fd04s stablerun.onM_fd04s];

rgcindices = 3;
[rasterrun conerun stablerun crs crsx ps resnorms residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', rgcindices, 'printpath', crprintpath, 'scaled_up', 10);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d07 d.d04s d.d08cm] = deal(cleared{:});


%% Unscaled C/R for Neuron reviewer 2
subplot(1,4,3)
plot(crsx(1,:)./3.*200, crs(:,:,3)', '.-', 'MarkerSize', 15);
title('ON Midget')
xlabel('% contrast')


%% Searching data007cm
if ~isfield(d, 'd07cm'), d.d07cm = load_data([piece '/data004_data007_data008-norefit/data007-from-data004_data007_data008/data007-from-data004_data007_data008'], loadopts); end
if ~isfield(d, 'd08cm'), d.d08cm = load_data([piece '/data004_data007_data008-norefit/data008-from-data004_data007_data008/data008-from-data004_data007_data008'], loadopts); end
if ~isfield(d, 'd04_07_08'), d.d04_07_08 = load_data([piece '/data004_data007_data008-norefit/data004_data007_data008-norefit'], loadopts); end

d.d07cm = read_stim_lisp_output(d.d07cm);
d.d07cm.stimulus = parse_stim_rgbs(d.d07cm.stimulus);
urgb = d.d07cm.stimulus.urgb;

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};

conerun = d.d04s;
rasterrun = d.d07cm;
map = map_ei(conerun, rasterrun);

% ON Parasols
celltype = 1;
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% There are some responses but cells look contaminated

% OFF Parasols
celltype = 2;
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% IDs from d04s
% Somewhat isolated: 6391
% Not well isolated: 602, 6619, 7201
% Weak: 2448, 4966, 4906
% Very weak: 3422, 4291, 5446
% Multiple surround 2671
% Weak surround: 5266
% Very weak surround: 3362

% ON Midgets
celltype = 3;
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = map{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end

%% data010 Allcones analysis
allconesprintpath = []; % printpath('allcones');

if ~isfield(d, 'd10'), d.d10 = load_data([piece '/data010'], loadopts); end
d.d10 = read_stim_lisp_output(d.d10, ':2011-12-13-2_f08_allcones');
d.d10.stimulus = parse_stim_rgbs(d.d10.stimulus);
d.d10.stimulus.triggers = d.d10.triggers(1:2:end);
cell_list_map = map_ei(d.d08s, d.d10, 'master_cell_type', {4});
d.d10.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, d.d08s.offM));
clear cell_list_map;

if ~isfield(d, 'd14s'), d.d14s = load_data([piece '/streamed/data014-0/data014-0'], staopts); end
cell_list_map = map_ei(d.d08s, d.d14s, 'master_cell_type', {4});
d.d14s.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, d.d08s.offM));
clear cell_list_map;

conerun = d.d08s;
datarun = d.d10;
stablerun = d.d14s;

conerun.rgcs = conerun.offM;
datarun.rgcs = datarun.offM_fd08s;
stablerun.rgcs = stablerun.offM_fd08s;

[ax rfweights crweights f gof p crs crsx] = allcones_plot_old(datarun, conerun, 'urgbs', d.d10.stimulus.urgbs(4), 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});
clear conerun datarun allconesprintpath


%% C/R data020cm
crprintpath = [];%printpath('cr', piece);

d.d20cm = load_data([piece '/data001-018-020-023-norefit/data020-from-data001_data018_data020_data023/data020-from-data001_data018_data020_data023'], loadopts);
d.d23cm = load_data([piece '/data001-018-020-023-norefit/data023-from-data001_data018_data020_data023/data023-from-data001_data018_data020_data023'], staopts);

% Map from d18s to concatenated (valid also for individuals from concatenated since norefit)
d.d01_18_20_23 = load_data([piece '/data001-018-020-023-norefit/data001-018-020-023-norefit'], loadopts);
cell_list_map = map_ei(d.d18s, d.d01_18_20_23);
mapped_offM = cell_list_map(get_cell_indices(d.d18s, d.d18s.offM));
mapped_onM  = cell_list_map(get_cell_indices(d.d18s, d.d18s.onM));
mapped_offM{3} = 1321;
mapped_offM{6} = 3586;
mapped_offM{7} = 3737;
mapped_offM{10} = 4862;
mapped_offM{11} = 5491;
mapped_offM{12} = 6901;
d.d01_18_20_23.rgcs = [mapped_offM mapped_onM];
clear cell_list_map mapped_offM mapped_onM;

conerun = d.d18s;
conerun.rgcs = [conerun.offM conerun.onM];

rasterrun = d.d20cm;
rasterrun.rgcs = d.d01_18_20_23.rgcs;

stablerun = d.d23cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun crs crsx ps resnorms residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 7, 'printpath', crprintpath, 'scaled_up', 10, 'YLim', [0 60]);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d20cm d.d18s d.d23cm] = deal(cleared{:});

% id3737, index 7, for GRI (1cone paper used different RF to have 2x2 for
% every RF in the paper)
x = 98.5; y = 153.5; del = 19;
axis([x x+del y y+del]);

% For J Freeman subunits analysis
save_crraster_data(conerun, rasterrun, fullfile(server_path, 'freeman', 'subunits_raster_data', piece), 'stablerun', stablerun, 'stablergcs', stablerun.rgcs);


%% Unscaled C/R for Neuron reviewer 2
subplot(1,4,1)
plot(crsx(1,:)./3.*200, crs(:,:,7)', '.-', 'MarkerSize', 15);
title('OFF Midget')
xlabel('% contrast')
ylabel('# spikes');


%% To get a 2x2 RF instead of 3x3 for the 1cone paper...
d.d11s = load_data(fullfile(piece, 'streamed/data011-0/data011-0'), staopts);
d.d20cm = load_data([piece '/data001-018-020-023-norefit/data020-from-data001_data018_data020_data023/data020-from-data001_data018_data020_data023'], loadopts);
d.d20cm = read_stim_lisp_output(d.d20cm);
datarun = d.d11s;
datarun = setsimplemarks(datarun, 3736, 3.5);
plot_rf_stimmap(d.d11s, 3736, d.d20cm.stimulus.mapnycpoly{1}, 'fit', false, 'scaled_up', 10);
axis([146 146+30 229 229+30]);

%% C/R data012cm/data013cm
crprintpath = [];%printpath('cr', piece);

if ~isfield(d, 'd12cm'), d.d12cm = load_data([piece '/data008_data012_data013_data014-norefit/data012-from-data008_data012_data013_data014/data012-from-data008_data012_data013_data014'], loadopts); end
if ~isfield(d, 'd13cm'), d.d13cm = load_data([piece '/data008_data012_data013_data014-norefit/data013-from-data008_data012_data013_data014/data013-from-data008_data012_data013_data014'], loadopts); end
if ~isfield(d, 'd14cm'), d.d14cm = load_data([piece '/data008_data012_data013_data014-norefit/data014-from-data008_data012_data013_data014/data014-from-data008_data012_data013_data014'], loadopts); end

if ~isfield(d.d12cm, 'offM_fd08s')
    if ~isfield(d, 'd08_12_13_14'), d.d08_12_13_14 = load_data([piece '/data008_data012_data013_data014-norefit/data008_data012_data013_data014-norefit'], loadopts); end
    cell_list_map = map_ei(d.d08s, d.d08_12_13_14);
    d.d12cm.onM_fd08s = cell_list_map(get_cell_indices(d.d08s, d.d08s.onM));
    d.d12cm.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, d.d08s.offM));
    d.d13cm.onM_fd08s = cell_list_map(get_cell_indices(d.d08s, d.d08s.onM));
    d.d13cm.offM_fd08s = cell_list_map(get_cell_indices(d.d08s, d.d08s.offM));    
    clear cell_list_map
end

conerun = d.d08s;
conerun.rgcs = [conerun.offM conerun.onM];
rasterrun = d.d12cm;
rasterrun.rgcs = [rasterrun.offM_fd08s rasterrun.onM_fd08s];
stablerun = d.d14cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun crs crsx ps resnorms residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all');
% Nice ones in here, a few seem to have failed xscale fitting for unclear
% reasons.  Need to add in the points from data013; more negative contrasts
% run there.


% For J Freeman subunits analysis
save_crraster_data(conerun, rasterrun, fullfile(server_path, 'freeman', 'subunits_raster_data', piece), 'stablerun', stablerun, 'stablergcs', stablerun.rgcs);


%% Continuing from above, something wrong with stimulus synchronization...

% length(d.d13cm.triggers) / 2
% length(d.d13cm.stimulus.pulses)
% 18 trials missing from triggers (i.e. 36 triggers missing since 2 per)

% hist(diff(d.d13cm.triggers)); % Distances between triggers look right, so not missing from the middle; probably missing from start?

rasterrun = d.d13cm;
rasterrun.rgcs = [rasterrun.offM_fd08s rasterrun.onM_fd08s];
[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, [zeros(18,1); rasterrun.triggers(1:2:end)], 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d13cm d.d18s d.d23cm] = deal(cleared{:});
% Sweetness; adding the 18 dummy triggers to the front fixes it!  
% Could try to recover those 18 trials with from negative trigger times?