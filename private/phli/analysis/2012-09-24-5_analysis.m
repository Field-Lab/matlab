% Basics
piece = '2012-09-24-5';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%% Streamed runs and cell picks
d.d01s = load_data([piece '/data001'], staopts);
d.d01s.onP = [336 3923 5509 6651];
d.d01s.onM = [2617 2977];

d.d03s = load_data([piece '/data003'], staopts);
d.d03s.onP = [3454 3739 5209 6049 6499 7237];
d.d03s.onM = [367 1011 2614 2855 3064 3919 4597 4971 5543 7042 7127];
d.d03s.onMsearched = 3650;
d.d03s.onMsurround = 3933;
d.d03s.offPsearched = [2521 2971 3616 722 4801 1459 887 2941 4456 3841 3647 4172];
d.d03s.offMsearched = [3152, 2851, 7231, 6931, 856, 5206, 5071, 3646, 3631, 346, 77, 2596, 4591, 4336, 3677];
d.d03s.offMsurround = [3031, 2866, 2776, 2686, 2492, 1336, 1186, 7081, 7351, 7576, 7561, 6706, 3031, 2492, 1186, 452, 451, 332, 5732, 5536, 5522, 5326, 5266, 5056, 4846, 4742, 4576, 3811, 3826, 3721, 3676, 3226, 7698, 7681, 7547 7576, 7561, 138, 122, 48, 827, 556, 541, 1111, 1051, 466, 751, 7576, 1382, 1291, 1217, 961, 1741, 1546, 1336, 1276, 1186, 827, 556, 541, 452, 451, 332, 1488, 1353, 1231, 1081, 3226, 4742, 4186, 5056, 4846, 4742, 4576, 4516 4339, 4337, 4306, 4051, 4066, 4576, 4516, 4411, 4339, 4337, 4306, 4051, 3826, 3721, 3691, 3676];


%        OFF parasols surround: 2102 & 2416 (not isolated)
%        OFF parasol surround: 3902 & 3841 & 3647(have to check cone-by-cone)
%        OFF parasols surround: 1501 (well isolated, not strong)
%        OFF parasols surround: 7726 (well isolated but not that strong), 691 & 7277 & 4188 (not isolated)
%        OFF parasol surround: 4188 fairly well isolated
%        OFF parasol surround: 3902, also surround from ON parasol 3739


%% Interdigitation
d.d01s = get_sta_summaries(d.d01s, {4}, 'robust_std_method', 5);
d.d03s = get_sta_summaries(d.d03s, {4}, 'robust_std_method', 5);
plot_rf_coloring(d.d01s, {4});
plot_rf_coloring(d.d03s, {4});


%% Likelihood surface figure
d.d01s = load_cones(d.d01s, 'bayes');
d.d03s = load_cones(d.d03s, 'bayes');

% Combining blurred things a little too much for my taste; just use d01
% Combine?
% expnt = 0.25;
% subplot(131);
% imagesc(matrix_scaled_up((norm_image(d.d01s.cones.dll)-0.5).^expnt,3));
% axis equal square
% subplot(132);
% imagesc(matrix_scaled_up((norm_image(d.d03s.cones.dll)-0.5).^expnt,3));
% axis equal square
% subplot(133);
% imagesc(matrix_scaled_up((norm_image(d.d01s.cones.dll + d.d03s.cones.dll)-0.5).^expnt,3));
% axis equal square

expnt = 0.25;
imagesc(matrix_scaled_up((norm_image(d.d01s.cones.dll)-0.5).^expnt,3));
axis equal square

% 100 um scale bar, based on Greg's measurement of the rig A OLED total screen size
plot([2500 2500-8.8235*3*10], [2600 2600], 'w')



%% C/R data006cm
crprintpath = [];%printpath('cr', piece);

if ~isfield(d, 'd06cm'), d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts); end
if ~isfield(d, 'd07cm'), d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts); end

if ~isfield(d, 'd03_06_07'), d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts); end
cell_list_map_d03s = map_ei(d.d03s, d.d03_06_07);

conerun = d.d03s;
conerun.rgcs = conerun.onM;
%conerun.rgcs = cell2mat(cellfun(@(name)(conerun.(name)), {'onP' 'onM' 'onMsearched'}, 'UniformOutput', false));
%conerun.rgcs = conerun.offPsearched;
%conerun.rgcs = conerun.offMsearched;

rasterrun = d.d06cm;
rasterrun.rgcs = cell_list_map_d03s(get_cell_indices(conerun, conerun.rgcs));

stablerun = d.d07cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun crs crsx ps resnorms residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d06cm d.d03s d.d07cm] = deal(cleared{:});
% ON midgets
%   2856, 3063, 4597 (d03s) beautiful, 97, 1010 (d03s) not bad, (but all low firing rates...)
% 
% Didn't use many negative contrasts, so most of the OFFs are not really usable :(.
% OFF parasols
%   722, 2971, 4801 (d03s) not bad
%
% OFF midgets
%   346, 2596, 2851, 3677, 4591, 6931 (d03s) not bad
%   3646, 5071 (d03s) soso
%

conerun.rgcs = conerun.offMsurround(1:10);

rasterrun = d.d06cm;
rasterrun.rgcs = cell_list_map_d03s(get_cell_indices(conerun, conerun.rgcs));

stablerun = d.d07cm;
stablerun.rgcs = rasterrun.rgcs;

[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'surround', true, 'printpath', crprintpath);
% OFF midgets surround
%   393, 2492, 2776, 2866, 3241, 3692 not bad (d03s)

if ~isfield(d, 'd05s'), d.d05s = load_data([piece '/data005'], loadopts); end
cell_list_map_d05s = map_ei(d.d05s, d.d03_06_07);



%% U/D data006 onM onP
udprintpath = [];

d.d06 = load_data([piece '/data006'], loadopts);
d.d07 = load_data([piece '/data007'], loadopts);

cell_list_map = map_ei(d.d03s, d.d06, 'master_cell_type', 'all');
d.d06.onM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onM));
d.d06.onP_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onP));
clear cell_list_map;

d.d06 = read_stim_lisp_output(d.d06);
d.d06.stimulus = parse_stim_rgbs(d.d06.stimulus);

conerun = d.d03s;
conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
conerun = set_polarities(conerun);
conerun = get_sta_fits_from_vision(conerun);
conerun = load_cones(conerun, 'acquisition');
conerun = make_mosaic_struct(conerun);
conerun = make_voronoi_masks(conerun);
conerun.wwrgcs = conerun.onM;

% stablerun = d.d07;
% stablerun = load_sta(stablerun, struct('load_sta', [], 'guess_stimulus', false));
% stablerun = set_polarities(stablerun);
% stablerun = get_sta_fits_from_vision(stablerun);
% stablerun.wwrgcs = stablerun.onM_fd03s;

ud06.wwrun = d.d06;
ud06.wwrun.wwrgcs = ud06.wwrun.onM_fd03s;
ud06.conerun = conerun;
% ud06.stablerun = stablerun;
ww_compound_plot(ud06, 'all', 'triggers', ud06.wwrun.triggers(1:2:end), 'printpath', udprintpath);

clear conerun stablerun udprintpath


%% U/D data006cm onM search from d05s
% A few new ones but not good
% udprintpath = [];
% 
% d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
% d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
% %d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts);
% 
% d.d05s = load_data([piece '/data005'], staopts);
% d.d03_06_07.mapd05s = map_ei(d.d05s, d.d03_06_07);
% 
% d.d06cm = read_stim_lisp_output(d.d06cm);
% d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);
% 
% conerun = d.d05s;
% conerun.wwrgcs = conerun.cell_types{3}.cell_ids;
% 
% ud06.conerun = conerun;
% ud06.wwrun = d.d06cm;
% ud06.wwrun.wwrgcs = d.d03_06_07.mapd05s(get_cell_indices(conerun, conerun.wwrgcs));
% % ud06.stablerun = stablerun;
% ww_compound_plot(ud06, 'all', 'triggers', ud06.wwrun.triggers(1:2:end), 'printpath', udprintpath);
% 
% clear conerun stablerun udprintpath


%% U/D data006cm onM onP, not sure if d03-06-07 has allspikesei; double check; can't be fixed now though due to loss of raw data :(
% I think raw data from these runs is now found actually!
udprintpath = [];
d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts);

cell_list_map = map_ei(d.d03s, d.d03_06_07, 'master_cell_type', 'all');
d.d03_06_07.onM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onM));
d.d03_06_07.onP_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onP));
d.d06cm.onM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onM));
d.d06cm.onP_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onP));
d.d07cm.onM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onM));
d.d07cm.onP_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onP));
clear cell_list_map;

d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);

conerun = d.d03s;
conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
conerun = set_polarities(conerun);
conerun = get_sta_fits_from_vision(conerun);
conerun = load_cones(conerun, 'acquisition');
conerun = make_mosaic_struct(conerun);
conerun = make_voronoi_masks(conerun);
conerun.wwrgcs = conerun.onM;

stablerun = d.d07cm;
stablerun = load_sta(stablerun, struct('load_sta', [], 'guess_stimulus', false));
stablerun = set_polarities(stablerun, 'cell_specs', {{1 3 5 6 7 8 9 10 11 12 13 14 15 16} {2}});
stablerun = get_sta_fits_from_vision(stablerun);
stablerun.wwrgcs = stablerun.onM_fd03s;

ud06cm.wwrun = d.d06cm;
ud06cm.wwrun.wwrgcs = ud06cm.wwrun.onM_fd03s;
ud06cm.conerun = conerun;
ud06cm.stablerun = stablerun;
ww_compound_plot(ud06cm, 'all', 'triggers', ud06cm.wwrun.triggers(1:2:end), 'printpath', udprintpath);

ud06cm.conerun.wwrgcs = conerun.onP;
ud06cm.stablerun.wwrgcs = stablerun.onP_fd03s;
ud06cm.wwrun.wwrgcs = ud06cm.wwrun.onP_fd03s;
ww_compound_plot(ud06cm, 'all', 'triggers', ud06cm.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear conerun stablerun udprintpath


%% data006cm: Any response from ON parasol cells when we double up the cone stimuli??
d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts);

cell_list_map = map_ei(d.d03s, d.d03_06_07, 'master_cell_type', 'all');
d.d03_06_07.onM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onM));
d.d03_06_07.onP_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onP));
d.d06cm.onM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onM));
d.d06cm.onP_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onP));
%d.d07cm.onM_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onM));
%d.d07cm.onP_fd03s = cell_list_map(get_cell_indices(d.d03s, d.d03s.onP));
clear cell_list_map;

d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);
urgb = d.d06cm.stimulus.urgb;

% Just plot all the on stimuli concatenated versus all the offs
% ons  = any(urgb.incr);
% offs = ~ons & ~urgb.blanks;
% for cellspec = [d.d06cm.onP_fd03s{:}]
%     figure
%     urgb_raster_subplot(d.d06cm, cellspec, d.d06cm.triggers(1:2:end), {{find(ons)} {find(offs)}});
% end

% More nuanced plot...
d.d03s = load_sta(d.d03s, 'load_sta', []);
d.d03s = set_polarities(d.d03s);
d.d03s = get_sta_fits_from_vision(d.d03s);
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellindex = [2 3 6]
    figure
    
    cellid = d.d06cm.onP_fd03s{cellindex};    
    urgb_raster_subplot(d.d06cm, cellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    cellid = d.d03s.onP(cellindex);
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false);
    set(gcf, 'Position', [516         186        1188         420]);
end

% Add this interesting ON parasol hit perfectly by ON midget stim, found by searching
cell_list_map = map_ei(d.d03s, d.d03_06_07);
cellid = 2988;
cellnum = d.d03s.cell_nums(cellid);
mappedcellid = cell_list_map{cellnum};
figure
urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
sanesubplot(2, 5, {1:2, 4:5});
plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false);
set(gcf, 'Position', [516         186        1188         420]);

for cellindex = [1 2 3 4 5 7]
    figure
    
    cellid = d.d06cm.onM_fd03s{cellindex};    
    urgb_raster_subplot(d.d06cm, cellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    cellid = d.d03s.onM(cellindex);
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false);
    set(gcf, 'Position', [516         186        1188         420]);    
end

% Add this interesting ON midget with surround response
cellid = 3933;
cellnum = d.d03s.cell_nums(cellid);
mappedcellid = cell_list_map{cellnum};
figure
urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
sanesubplot(2, 5, {1:2, 4:5});
plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false);
set(gcf, 'Position', [516         186        1188         420]);


%% data006cm This is actually good for searching too...
% Didn't include ON midget 2614 in original search; should go back and look
% at parasols for this one.
d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts);

d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);
urgb = d.d06cm.stimulus.urgb;

cell_list_map_d03s = map_ei(d.d03s, d.d03_06_07);

d.d03s = load_sta(d.d03s, 'load_sta', []);
d.d03s = set_polarities(d.d03s);
d.d03s = get_sta_fits_from_vision(d.d03s);
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

% Search all ON Parasols
celltype = 1;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = d.d03s.cell_types{celltype}.cell_ids
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false);
end
% Interesting one found!     2988 (matching ON midget: 2855, matching OFF parasols: )
% Slightly interesting ones: 2768 3920

% Search all OFF Parasols
celltype = 2;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = d.d03s.cell_types{celltype}.cell_ids
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Lots of interesting, if hard to work with, stuff in here!!
%
% Lots of costimulated, good number of surround stimulated too!  But not
% many good cases of isolated stimulation, especially surround stimulation
%
% Especially nice ones, with overlap on ON midgets or even ON parasols:
%   ON parasol 2988 (only fires with double cone ON) 
%        ON midget 2855 is pretty, 
%        OFF parasols: 2521 (reasonably isolated), 2971 (not isolated)
%        OFF parasols surround: 2102 & 2416 (not isolated)
%   ON parasol 3739
%        OFF parasol: 3616
%        OFF parasol surround: 3902 & 3841 & 3647(have to check cone-by-cone)
%   ON parasol 7237
%        OFF parasol: 722 (not completely isolated)
%   ON parasol 5209
%        OFF parasol: 4801 (not completely isolated, check cone-by-cone)
%   ON midget 1011
%        OFF parasols: 1459 (reasonably isolated, edge), 887 (not completely isolated)
%        OFF parasols surround: 1501 (well isolated, not strong)
%   ON midget 367
%        OFF parasols surround: 7726 (well isolated but not that strong), 691 & 7277 & 4188 (not isolated)
%   ON midget 3064
%        OFF parasol edge: 2941 (not well isolated; interesting to look cone-by-cone)
%   ON midget 4597
%        OFF parasol: 4456 well isolated
%        OFF parasol surround: 4188 fairly well isolated
%   ON midget surround 3933
%        OFF parasol: 3841 & 3647, also surround from ON parasol 3739; check cone-by-cone
%        OFF parasol edge: 4172, also surround from ON midget 4597 and others
%        OFF parasol surround: 3902, also surround from ON parasol 3739

% Search all ON Midgets
celltype = 3;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = d.d03s.cell_types{celltype}.cell_ids
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Not too much here; mostly just the intended ones hit
%   ON parasol 3739
%       ON midget 3650 pretty decent for a few of the cones
%   ON midget 6802, 6197, very weak
%   ON midget 3336 not really any good
%
% ON midget 3933 shows surround stimulation so that's cool; should search back for the region in OFF parasols.

% Search all OFF Midgets
celltype = 4;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = d.d03s.cell_types{celltype}.cell_ids(end)
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end
% Hoo boy...  175 mapped OFF midgets (201 in the streamed data).  Lots of
% crazy stuff in here.  Definitely need to check spike sorting...
%
%   ON parasol 2988 & ON midget 2855
%       OFF midget edge: 3152, 2851, 2641
%       OFF midget surrounds: 2401
%       OFF midget surrounds shared: 3031, 2866, 2776, 2686, 2492, 2341, 2251, 1741, 1546, 1336, 1186 %%
%   ON parasol 7237
%       OFF midget center: 7231, 6931
%       OFF midget edge: 856
%       OFF midget surrounds (may have mixed up some shared in here): 7081, 7351, 7336, 7576, 7561, 6706
%       OFF midget surrounds shared: 3031, 2492, 1186, 452, 451, 332
%   ON parasol 5209
%       OFF midget center: 5206
%       OFF midget edge: 5191, 5071
%       OFF midget surrounds: 5732, 5536, 5522, 5326, 5266
%       OFF midget surrounds shared: 5056, 4846, 4742, 4576
%   ON parasol 3739
%       OFF midget center+surround shared: 3646, 3631
%       OFF midget surrounds: 3811
%       OFF midget surrounds shared: 3826, 3721, 3676, 3226
%   ON midget 367
%       OFF midget center: 346 (perfect for U/D)
%       OFF midget edge: 421, 77
%       OFF midget surrounds: 7698, 7681, 7547 (well isolated), 7576, 7561, 138, 122, 48 %%
%       OFF midget surrounds shared: 827, 556, 541
%       OFF midget, not much despite being closer: 7756
%       OFF midget, not much despite being seemingly close enough: 7516, 751 %%
%   ON midget 1011
%       OFF midget edge: 1111, 1051, 466
%       OFF midget edge shared: 751
%       OFF midget surrounds: 7576, 1382, 1291, 1217, 961
%       OFF midget surrounds shared: 1741, 1546, 1336, 1276, 1186, 827, 556, 541, 452, 451, 332 %%
%       OFF midget weak surround: 1488, 1353, 1231, 1081
%   ON midget 2614
%       OFF midget center: 2596 (maybe usable for OFF midget U/D)
%   ON midget 3064
%       OFF midget edge: 3211
%       OFF midget surrounds shared: 3226
%   ON midget 4597
%       OFF midget center: 4591
%       OFF midget edge: 4336
%       OFF midget surrounds: 4742, 4186
%       OFF midget surrounds shared: 5056, 4846, 4742, 4576, 4516 4339, 4337, 4306, 4051 %%
%   ON midget surround 3933
%       OFF midget center: 3677 (perfect; can use for OFF U/D)
%       OFF midget edge: 4141, 4036
%       OFF midget surround: 4066
%       OFF midget surrounds shared: 4576, 4516, 4411, 4339, 4337, 4306, 4051, 3826, 3721, 3691, 3676 %%
%   Nulls, i.e. perhaps too far for surround stim (but could also just be shitty cell...)
%       OFF midget 1951 (were couple more just after this), 1786, 1771, 1756 %%

% SBCs
celltype = 5;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = d.d03s.cell_types{celltype}.cell_ids
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end


%% data006cm Revisit search using mapping from d05s; in particular, many more ON midgets there
d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts);

d.d05s = load_data([piece '/data005'], loadopts);
cell_list_map_d05s = map_ei(d.d05s, d.d03_06_07);

d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);
urgb = d.d06cm.stimulus.urgb;

d.d05s = load_sta(d.d05s, 'load_sta', []);
d.d05s = set_polarities(d.d05s);
d.d05s = get_sta_fits_from_vision(d.d05s);

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d05s;
rasterrun = d.d06cm;

% Not much for ON Parasols

% Search all OFF Parasols
celltype = 2;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Seems like some new ones in here:
%   ON parasol 2988 & ON midget 2855 (d03s)
%       OFF parasol shared surround 1907 (d05s)
%
%   ON parasol 3739 (d03s)
%       OFF parasol edge 3616 (d05s)
%       OFF parasol shared edge 3857 (d05s)
%       OFF parasol shared surround weak 3901 (d05s)
%
%   ON parasol 5209 (d03s super weak)
%       OFF parasol shared edge 4666 (d05s)
%
%   ON parasol 7237 (d03s)
%       OFF parasol edge 6721, 7381 (d05s)
%       OFF parasol shared surround 7126 (d05s)
%
%   ON midget 367 (d03s)
%       OFF parasol center 122 (d05s)
%       OFF parasol edge 422 (d05s)
%       OFF parasol shared surround 766, 7126 (d05s)
%       OFF parasol isolated surround 7726 (d05s)
%
%   ON midget 1011 (d03s)
%       OFF parasol edge 842, 1277 (d05s)
%       OFF parasol shared surround 766, 1907 (d05s)
%
%   ON midget 2614 (d03s)
%       OFF parasol shared surround 1907 (d05s)
%
%   ON midget 3933 (d03s)
%       OFF parasol edge 4261 (d05s)
%       OFF parasol shared edge 3857 (d05s)
%       OFF parasol shared surround weak 3901 (d05s)
%
%   ON midget 4597 (d03s)
%       OFF parasol center 4456 (d05s)
%       OFF parasol shared edge 4666 (d05s)
%
%   OFF parasol center 5146, 5941, 6091 (d05s)
%   OFF parasol edge 4352, 6574 (d05s)


% Search all ON Midgets
celltype = 3;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Couple new ones
%   ON parasol 3739 (d03s)
%       ON midget edge 3738 (d05s)
%       ON midget shared surround 3706 (d05s)
%
%   ON midget 3933 (d03s)
%       ON midget shared surround 3706 (d05s)
%
% ON midget edge 6346 (d05s)
% ON midget edge weak 6803 (d05s)
% ON midget edge 3650 (d03s; already found) -> 3648 d05s
% ON midget shared surround 3933 (d03s; already found) -> 3931 d05s
%


%% Condensed rasters for specific interesting groups of cells
% (pick set here and then run next cell to plot)

% All four types (weak ON parasol)
ons = [2988 2855];
ons_d05s = [];
offs = [2521 2416 2102 2851 2341];
offs_d05s = [];
%offs = [2521 2971 2102 2416]; % P
%offs = [offs 3152, 2851, 2641, 2401]; % M
%offs = [offs 3031, 2866, 2776, 2686, 2492, 2341, 2251, 1741, 1546, 1336, 1186]; % M shared

% All four types (weak ON parasol)
ons = [3739 3650];
ons_d05s = [3738 3706];
offs = [3616 3631     3646 3902 3826 3811];
offs_d05s = [];
%offs = [3616 3902 3841 3647   3646 3631 3811 3826 3721 3676 3226];

% Divergence
% ons = []; ons_d05s = []; offs = []; offs_d05s = [];
ons = [4597]; ons_d05s = []; offs = [4591]; offs_d05s = []; % ON and OFF M, only one strong cone in OFF, but looks pretty good.
ons = []; ons_d05s = []; offs = [3677 4172]; offs_d05s = []; % OFF M and P, pretty decent
% ons = [367]; ons_d05s = []; offs = [346 77]; offs_d05s = []; % ON and OFF M, but ON M only ~3 Hz.
% ons = [1011]; ons_d05s = []; offs = [1111 1051 466 751]; offs_d05s = [842 1277];
% ons = [2614]; ons_d05s = []; offs = [2596]; offs_d05s = [1907];
% ons = [3064]; ons_d05s = []; offs = [3211 2941]; offs_d05s = [];


% Surprising surround effects
ons = [367];
ons_d05s = [];
offs = [7726 691 7277          346 77         48 421 122 138   7756   7561 7681 7698 7547 7516];
offs_d05s = [];
%offs = [7726 691 7277 4188            346 421, 77 7698, 7681, 7547, 7576, 7561, 138, 122, 48  827, 556, 541  7756 7516 751];

% Nice surround gradient
ons = [1011];
ons_d05s = [];
offs = [1051   1111, 466    1081, 1276, 1217, 961, 1353, 1291, 1382, 1546, 1741, 1231, 1488];
offs_d05s = [];
%offs = [1111, 1051, 466       751         7576, 1382, 1291, 1217, 961       1741, 1546, 1336, 1276, 1186, 827, 556, 541, 452, 451, 332     1488, 1353, 1231, 1081];
% Look back through parasols

% ON midget surround
ons = [3933];
ons_d05s = [];
offs = []; % FILL IN
offs_d05s = [];


%% Condensed rasters for specific interesting groups of cells
% (pick set above and then run this to plot)
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = ons;
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2, 'hist_line_width', 2, 'hist_color', 'm');
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
    set(gcf, 'Position', [516         186        1188         420]);
end
for cellid = ons_d05s;
    cellnum = d.d05s.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2, 'hist_line_width', 2, 'hist_color', 'm');
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d05s, cellid, d.d06cm.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
    set(gcf, 'Position', [516         186        1188         420]);
end
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = offs;
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2, 'hist_line_width', 2, 'hist_color', 'm');
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
    set(gcf, 'Position', [516         186        1188         420]);
end
for cellid = offs_d05s;
    cellnum = d.d05s.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2, 'hist_line_width', 2, 'hist_color', 'm');
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d05s, cellid, d.d06cm.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
    set(gcf, 'Position', [516         186        1188         420]);
end


%% data004cm; do some searching here too
d.d01_04_05 = load_data([piece '/d01-04-05-norefit/d01-04-05-norefit'], loadopts);
d.d04cm = load_data([piece '/d01-04-05-norefit/data004/data004'], loadopts);
d.d05cm = load_data([piece '/d01-04-05-norefit/data005/data005'], loadopts);

cell_list_map_d01s = map_ei(d.d01s, d.d01_04_05);

d.d04cm = read_stim_lisp_output(d.d04cm);
d.d04cm.stimulus = parse_stim_rgbs(d.d04cm.stimulus);
urgb = d.d04cm.stimulus.urgb;

d.d01s = load_sta(d.d01s, 'load_sta', []);
d.d01s = set_polarities(d.d01s);
d.d01s = get_sta_fits_from_vision(d.d01s);
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d01s;
rasterrun = d.d04cm;


% Search all ON Parasols
celltype = 1;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false);
end
% ON parasol 336 responds to double cones


% Search all OFF Parasols
celltype = 2;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = [7098, 6106, 4246, 3902]%conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Multiple region OFF parasol surrounds
%   7278, 5990, 5853, 5747, 4801, 3618, 2822, 2748, 2042, 1324
%
% Cusp OFF parasol surround
%   1501, 1381, 46 (ON parasol 336 region)
%
% One region dominant OFF parasol surround
%   7098, 6106, 4246, 3902
%
% Edge OFF parasol
%   5043, 4172, 3724, 2268, 1188 & 261 (ON parasol 336 region)
%
% Mostly isolated OFF parasol center
%   3826, 2521, 887 (inside cusp?)
%
% Non-isolated OFF parasol center
%   6886, 6858, 6423
%
% Not much stimulated OFF parasol
%   6797, 5896, 4591, 692, 151


% Search all ON Midgets
celltype = 3;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Edge: 3933
% Center: 2977
% Should search for more


% Search all OFF Midgets
celltype = 4;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end
% Center OFF midgets
%   466 (ON parasol 336 region)
%   6558, 3781, 3287, 2852, 2596
%
% Cusp OFF midgets
%   2584
%
% Edge OFF midgets
%   335, 332 (ON parasol 336 region)
%   6841, 6664, 6631, 6601, 6511, 6496, 6392, 6391, 6256, 5626, 3676, 2642, 2641, 2431, 2192, 871
%
% Isolated surround OFF midgets:
%   7711, 828, 529, 347, 138, 122 (ON parasol 336 region)
%   707, 423 (ON parasol 336 region weak)
%   6466, 4173, 4127, 4112, 4051, 3056, 3828, 3811, 2776, 2747
%   3797, 3736 (very weak)
%
% Mostly isolated surround OFF midgets:
%   7681 (ON parasol 336 region)
%
% Dominated surround OFF midgets
%   7578, 7561, 7366, 1111, 1081, 1051, 916, 886, 796, 451 (ON parasol 336 region)
%   7396, 7264, 7248, 7247, 7096, 7081, 6991, 6856, 6856, 6811, 6707, 6560, 6272, 5986, 5716, 5582, 5506, 5492, 5237, 5209, 4771, 4411, 4338, 4143, 4142, 3721, 3646, 3632, 3316, 2356, 1741, 1727, 1338 %
%   6946, 6093 (low baseline)
%   5297, 4592, 4186, 1982 (weak)
%   5267, 5073, 5071, 4726, 4516, 4340, 4336, 2028, 1892, 1681, 1231i, 1173 (very weak)
%
% Multiple surround OFF midgets
%   7548, 7471, 7187, 7173, 7158, 7141, 6931, 6917, 6587, 6303, 6302, 6287, 6286, 5611, 4966, 4876, 4711, 4501, 4486, 4233, 4216, 4007, 3692, 3466, 3455, 3406, 3151, 3121, 3106, 3046, 3031, 2866, 2686, 2566, 2492, 2401, 2252, 2191, 1471, 856 %
%   7516, 7501, 6031, 5943, 5881, 5807, 5672, 5416, 4847, 4846, 3616, 3362, 3226, 2386, 1801, 1666, 1531, 1516, 1277 (weak) (??? possibly mislabeled)
%
% Not much stimulated
%   5671, 5657, 5401, 5329, 5327, 5191, 5056, 4398, 3061, 1891, 1771, 1756, 1637, 1486, 1426, 1367
%   2056, 1307 (fairly close)
%   961, 49, 34 (fairly close ON parasol 336 region)
%
% Weird OFF midgets
%   7232, 6062, 5522, 5178


% Search all SBCs
celltype = 5;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end

    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);

    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% 2313 shows some nice surround.