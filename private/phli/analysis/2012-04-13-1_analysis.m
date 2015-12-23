%% Basics
piece = '2012-04-13-1';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
udprintpath = [];
keep_vars = {'piece'; 'loadopts'; 'udprintpath'; 'allconesprintpath'; 'keep_vars'; 'd'};


%% The runs with picked cells
streamedopts = loadopts;
streamedopts.load_sta = true;
streamedopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
streamedopts.set_polarities = {'guess', true};

d.d02s = load_data([piece '/streamed/data002/data002'], streamedopts);
d.d02s.offM = [288 1321 2118 2389 2941 3001 5866 6286 6395 6663 6857 7700];

d.d06s = load_data([piece '/streamed/data006/data006'], streamedopts);
d.d06s.offM = [2536 5927 6136 7700];


%% Mapping individual
d.d02 = load_data([piece '/data002'], loadopts);
cell_list_map = map_ei(d.d02s, d.d02, 'master_cell_type', {4});
d.d02.offM_fd02s = cell_list_map(get_cell_indices(d.d02s, [d.d02s.offM]));
clear cell_list_map;
% Looks good except there are duplicates for 5867: 5868 5869
% 426 tiny tiny bit contamination?
% 1683 tiny bit of contamination?
% 2941 not very clean
% 5867 split as expected
% 6664, some contamination
% 6751 not so clean
% 7700 okay, little contamination, maybe lost spikes to 7699 7701 7703

d.d02rch = load_data([piece '/data002-from-data002CH'], loadopts);
cell_list_map = map_ei(d.d02s, d.d02rch, 'master_cell_type', {4});
d.d02rch.offM_fd02s = cell_list_map(get_cell_indices(d.d02s, [d.d02s.offM]));
clear cell_list_map;
% 1683 tiny bit of contamination?
% 2941 not clean
% 6664, some contamination but looks mostly okay
% 6751 not so clean
% 7700 okay, little contamination, maybe lost spikes to 7699 7703 7704

d.d04 = load_data([piece '/data004'], loadopts);
cell_list_map = map_ei(d.d02s, d.d04, 'master_cell_type', 'all');
d.d04.offM_fd02s = cell_list_map(get_cell_indices(d.d02s, [d.d02s.offM]));
clear cell_list_map;

d.d04rch = load_data([piece '/data004-from-data002CH'], loadopts);
cell_list_map = map_ei(d.d02s, d.d04rch, 'master_cell_type', 'all');
d.d04rch.offM_fd02s = cell_list_map(get_cell_indices(d.d02s, [d.d02s.offM]));
clear cell_list_map;

d.d05 = load_data([piece '/data005'], loadopts);
cell_list_map = map_ei(d.d02s, d.d05, 'master_cell_type', 'all');
d.d05.offM_fd02s = cell_list_map(get_cell_indices(d.d02s, [d.d02s.offM]));
clear cell_list_map;

d.d05rch = load_data([piece '/data005-from-data002CH'], loadopts);
cell_list_map = map_ei(d.d02s, d.d05rch, 'master_cell_type', 'all');
d.d05rch.offM_fd02s = cell_list_map(get_cell_indices(d.d02s, [d.d02s.offM]));
clear cell_list_map;

d.d06 = load_data([piece '/data006'], loadopts);
cell_list_map = map_ei(d.d02s, d.d06, 'master_cell_type', {4});
d.d06.offM_fd02s = cell_list_map(get_cell_indices(d.d02s, [d.d02s.offM]));
clear cell_list_map;

d.d09 = load_data([piece '/data009'], loadopts);
cell_list_map = map_ei(d.d06s, d.d09, 'master_cell_type', {4});
d.d09.offM_fd06s = cell_list_map(get_cell_indices(d.d06s, [d.d06s.offM]));
clear cell_list_map;

d.d10 = load_data([piece '/data010'], loadopts);
cell_list_map = map_ei(d.d06s, d.d10, 'master_cell_type', {4});
d.d10.offM_fd06s = cell_list_map(get_cell_indices(d.d06s, [d.d06s.offM]));
clear cell_list_map;


%% C/R d04
crprintpath = [];%printpath('cr', piece);

conerun = d.d02s;
conerun.rgcs = conerun.offM;
rasterrun = d.d04;
rasterrun.rgcs = rasterrun.offM_fd02s;
stablerun = d.d06;
stablerun.rgcs = stablerun.offM_fd02s;

[rasterrun conerun stablerun] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rgcindices', 'all', 'printpath', crprintpath);
cleared = cellfun(@(s)(rmfield(s, 'rgcs')), {rasterrun conerun stablerun}, 'UniformOutput', false);
[d.d04 d.d02s d.d06] = deal(cleared{:});


%% U/D d05, new for subunits paper

d.d04 = read_stim_lisp_output(d.d04);
d.d05 = read_stim_lisp_output(d.d05);
d.d04.stimulus = parse_stim_rgbs(d.d04.stimulus);
d.d05.stimulus = parse_stim_rgbs(d.d05.stimulus);

conerun = d.d02s;
conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
conerun = set_polarities(conerun);
conerun = get_sta_fits_from_vision(conerun);
conerun = load_cones(conerun, 'acquisition');
conerun = make_mosaic_struct(conerun);
conerun = make_voronoi_masks(conerun);
conerun.wwrgcs = conerun.offM;

stablerun = d.d06;
stablerun = load_sta(stablerun, struct('load_sta', [], 'guess_stimulus', false));
stablerun = set_polarities(stablerun);
stablerun = get_sta_fits_from_vision(stablerun);
stablerun.wwrgcs = stablerun.offM_fd02s;

rasterrun = d.d05;
figure
outstructs = ud_compound_plot(conerun, rasterrun, 2389, 'triggers', rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rfopts', {'scaled_up', 10}, 'rasteropts', {'axopts', {{'Box', 'off', 'XTick', [0 0.6]}}});
outstructs = cell2mat(invertselect(outstructs, @isempty));
AX = vertcat(outstructs.AX);
set(AX(:,1), 'YTick', []);
set(AX(:,2), 'YTick', [0 40]);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(:,1), 'XTickLabel', []);
set(AX(setdiff(1:20, 4:4:20),2), 'XTickLabel', []);
set(gcf, 'Position', [2148 69 1185 846]);

figure
outstructs = ud_compound_plot(conerun, rasterrun, 3001, 'triggers', rasterrun.triggers(1:2:end), 'stabilityrun', stablerun, 'rfopts', {'scaled_up', 10}, 'rasteropts', {'axopts', {{'Box', 'off', 'XTick', [0 0.6]}}});
outstructs = cell2mat(invertselect(outstructs, @isempty));
AX = vertcat(outstructs.AX);
set(AX(:,1), 'YTick', []);
set(AX(:,2), 'YTick', [0 40]);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(1:end-3,2), 'YTickLabel', []);
set(AX(:,1), 'XTickLabel', []);
set(AX(setdiff(1:20, 4:4:20),2), 'XTickLabel', []);
set(gcf, 'Position', [2148 69 1185 846]);


%% U/D d04 d05 d04rch d05rch

d.d04 = read_stim_lisp_output(d.d04);
d.d05 = read_stim_lisp_output(d.d05);
d.d04.stimulus = parse_stim_rgbs(d.d04.stimulus);
d.d05.stimulus = parse_stim_rgbs(d.d05.stimulus);

%ww_raster_plot(d.d04, s.s04, offM_fd02s_d04, -1, 'triggers', d.d04.triggers(1:2:end))
%ww_raster_plot(d.d05, s.s05, offM_fd02s_d05, -1, 'triggers', d.d05.triggers(1:2:end))

conerun = d.d02s;
conerun = load_sta(conerun, struct('load_sta', [], 'guess_stimulus', false));
conerun = set_polarities(conerun);
conerun = get_sta_fits_from_vision(conerun);
conerun = load_cones(conerun, 'acquisition');
conerun = make_mosaic_struct(conerun);
conerun = make_voronoi_masks(conerun);
conerun.wwrgcs = conerun.offM;

stablerun = d.d06;
stablerun = load_sta(stablerun, struct('load_sta', [], 'guess_stimulus', false));
stablerun = set_polarities(stablerun);
stablerun = get_sta_fits_from_vision(stablerun);
stablerun.wwrgcs = stablerun.offM_fd02s;

ud04.wwrun = d.d04;
ud04.wwrun.wwrgcs = ud04.wwrun.offM_fd02s;
ud04.conerun = conerun;
ud04.stablerun = stablerun;
ww_compound_plot(ud04, 'all', 'triggers', ud04.wwrun.triggers(1:2:end), 'printpath', udprintpath);
clear conerun stablerun

% For J Freeman subunits analysis
save_udraster_data(ud04, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


d.d04rch = read_stim_lisp_output(d.d04rch);
d.d04rch.stimulus = parse_stim_rgbs(d.d04rch.stimulus);
ud04rch = ud04;
ud04rch.wwrun = d.d04rch;
ud04rch.wwrun.wwrgcs = ud04rch.wwrun.offM_fd02s;
ww_compound_plot(ud04rch, 'all', 'triggers', ud04rch.wwrun.triggers(1:2:end), 'printpath', udprintpath);


% For J Freeman subunits analysis
save_udraster_data(ud04rch, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


ud05 = ud04;
ud05.wwrun = d.d05;
ud05.wwrun.wwrgcs = ud05.wwrun.offM_fd02s;
ww_compound_plot(ud05, 'all', 'triggers', ud05.wwrun.triggers(1:2:end), 'printpath', udprintpath);


% For J Freeman subunits analysis
save_udraster_data(ud05, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


d.d05rch = read_stim_lisp_output(d.d05rch);
d.d05rch.stimulus = parse_stim_rgbs(d.d05rch.stimulus);
ud05rch = ud05;
ud05rch.wwrun = d.d05rch;
ud05rch.wwrun.wwrgcs = ud05rch.wwrun.offM_fd02s;
ww_compound_plot(ud05rch, 'all', 'triggers', ud05.wwrun.triggers(1:2:end), 'printpath', udprintpath);


% For J Freeman subunits analysis
save_udraster_data(ud05rch, fullfile(server_path, 'freeman', 'subunits_raster_data', piece));


clear printpath


%% Separate for figures
ww_compound_plot(ud04rch, 4,  'triggers', ud04rch.wwrun.triggers(1:2:end), 'set_ylims', [0 34], 'scaled_up', 10, 'MarkerSize', 15, 'FontSize', 14);
ww_compound_plot(ud05,    10, 'triggers', ud05.wwrun.triggers(1:2:end),    'set_ylims', [0 16], 'scaled_up', 10, 'MarkerSize', 15, 'FontSize', 14);
ww_compound_plot(ud05rch, 4,  'triggers', ud05rch.wwrun.triggers(1:2:end), 'set_ylims', [0 22], 'scaled_up', 10, 'MarkerSize', 15, 'FontSize', 14);


%% Mapping concatenated
d.d02cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data002-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data002-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);

cell_list_map = map_ei(d.d02s, d.d02cm, 'master_cell_type', 'all');
d.d02cm.offM_fd02s = cell_list_map(get_cell_indices(d.d02s, [d.d02s.offM]));
clear cell_list_map;
d.d02cm.offM_fd02s{1} = 286;
d.d02cm.offM_fd02s{10} = 6799;

cell_list_map = map_ei(d.d06s, d.d02cm, 'master_cell_type', 'all');
d.d02cm.offM_fd06s = cell_list_map(get_cell_indices(d.d06s, [d.d06s.offM]));
clear cell_list_map;
d.d02cm.offM_fd06s{3} = 6286;


%% Timecourse templates
% d.d00cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data000-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data000-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
% d.d01cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data001-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data001-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
% d.d02cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data002-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data002-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
% dnames =  = {'d00cm' 'd01cm' 'd02cm'};

d.d01 = load_data([piece '/data001'], loadopts);
dnames = {'d01'};

for dname = dnames
    datarun = d.(dname{1});
    datarun = load_sta(datarun, 'load_sta', []);
    
    for cell_type = [1 2 3 4]
        template = zeros(datarun.stas.depth,1);
        cnums = get_cell_indices(datarun, {cell_type});
        for cnum = cnums
            timecourse = datarun.vision.timecourses(cnum);
            template = template + timecourse.r;
            template = template + timecourse.g;
            template = template + timecourse.b;
        end; clear cnum cnums timecourse;
        datarun.vision.timecourse_templates{cell_type} = template;
    end
    
    d.(dname{1}) = datarun; clear datarun;
end
leave(keep_vars{:});


%% Timecourse C/R d04 d05
d.d04 = read_stim_lisp_output(d.d04);
d.d04.stimulus = parse_stim_rgbs(d.d04.stimulus);
d.d05 = read_stim_lisp_output(d.d05);
d.d05.stimulus = parse_stim_rgbs(d.d05.stimulus);

d.d04rch = read_stim_lisp_output(d.d04rch);
d.d04rch.stimulus = parse_stim_rgbs(d.d04rch.stimulus);
d.d05rch = read_stim_lisp_output(d.d05rch);
d.d05rch.stimulus = parse_stim_rgbs(d.d05rch.stimulus);

% d.d04cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data004-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data004-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
% d.d05cm = load_data([piece
% '/d00-01-02-04-05-06-07-09-10-norefit/data005-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data005-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);

trname = 'd01'; % Template run
template_index = 4;
drname = 'd04rch';
rgcsfield = 'offM_fd02s';
template_fit_1cone_flash_contrast_response_calc(d.(drname), rgcsfield, d.(trname), template_index)

leave(keep_vars{:});


%% C/R simulation

d.d04 = read_stim_lisp_output(d.d04);
d.d04.stimulus = parse_stim_rgbs(d.d04.stimulus);
flashstim = d.d04.stimulus;

templaterun = d.d01;
template = [0 fliplr(templaterun.vision.timecourse_templates{4}')];
[flashtemplate, flashtemplatex] = flash_timecourse_template(template, templaterun, flashstim);

polarity = -1;
baseline = 1;
basic_sim_spikes(template.*polarity, 0.001, 50, 'baseline', baseline, 'xstart', -1);

leave(keep_vars{:})


%% data007 gratings analysis
if ~exist('s') || ~isfield(s, 's07')
    s.s07 = read_stim_lisp_output(piece, 's07');
end



%% data009 Allcones analysis
allconesprintpath = []; %[printpath('/allcones')];

if ~exist('d') || ~isfield(d, 'd09cm')
    d.d09cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data009-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data009-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
end
d.d09cm = read_stim_lisp_output(d.d09cm, ':2012-04-13-1:f06_allcones');
d.d09cm.stimulus = parse_stim_rgbs(d.d09cm.stimulus);

conerun = d.d06s;
datarun = d.d09cm;
conerun.rgcs = conerun.offM;
datarun.rgcs = d.d02cm.offM_fd06s;
datarun.stimulus.triggers = datarun.triggers(1:2:end);

conerun = load_sta(conerun, 'load_sta', []);
conerun = get_sta_fits_from_vision(conerun);
conerun = set_polarities(conerun);

if ~isfield(d, 'd10cm')
    d.d10cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data010-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data010-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
end
stablerun = d.d10cm;
stablerun.rgcs = d.d02cm.offM_fd06s;
stablerun = load_sta(stablerun, 'load_sta', []);

allcones_plot(datarun, conerun, 'printpath', allconesprintpath, 'urgbs', datarun.stimulus.urgbs(2:4), 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});


% Alternative
d.d09 = read_stim_lisp_output(d.d09, ':2012-04-13-1:f06_allcones');
d.d09.stimulus = parse_stim_rgbs(d.d09.stimulus);
datarun = d.d09;
datarun.rgcs = d.d09.offM_fd06s;
datarun.stimulus.triggers = datarun.triggers(1:2:end);
allcones_plot(datarun, conerun, 'printpath', allconesprintpath, 'urgbs', datarun.stimulus.urgbs(2:4), 'stablerun', stablerun, 'stimmasklocal', true, 'rfstimmapopts', {'scaled_up', 10});


leave(keep_vars{:});


%%


leave(keep_vars{:});


%% Mapping concatenated
d.d02s = load_data([piece '/streamed/data002/data002'], loadopts);
d.d02cm = load_data([piece '/d00-01-02-04-05-06-10-norefit/data002-from-data000_data001_data002_data004_data005_data006_data010/data002-from-data000_data001_data002_data004_data005_data006_data010'], loadopts);

cell_list_map = map_ei(d.d02s, d.d02cm, 'master_cell_type', 'all');
offM_fd02s_d02cm = cell_list_map(get_cell_indices(d.d02s, [offM_d02s]));
clear cell_list_map;
offM_fd02s_d02cm{1} = 286;
offM_fd02s_d02cm{10} = 6799;


%% U/D d04cm d05cm

if ~isfield(s, 's04')
    s.s04 = read_stim_lisp_output(piece, 's04');
    s.s04 = parse_stim_rgbs(s.s04);
end
if ~isfield(s, 's05')
    s.s05 = read_stim_lisp_output(piece, 's05');
    s.s05 = parse_stim_rgbs(s.s05);
end

d.d04cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data004-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data004-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
ww_raster_plot(d.d04cm, s.s04, offM_fd02s_d02cm, -1, 'triggers', d.d04cm.triggers(1:2:end))

d.d05cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data005-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data005-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
ww_raster_plot(d.d05cm, s.s05, offM_fd02s_d02cm, -1, 'triggers', d.d05cm.triggers(1:2:end))


%%
d.d07cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data007-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data007-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
s.s07 = read_stim_lisp_output(piece, 's07');
