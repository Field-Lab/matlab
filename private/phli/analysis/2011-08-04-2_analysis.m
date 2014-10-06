%% Load data
piece = '2011-08-04-2';

d01str = load_data('2011-08-04-2/streamed/data001/data001');
d05str = load_data('2011-08-04-2/streamed/data005-0/data005-0');
d06 = load_data('2011-08-04-2/data006');
d07 = load_data('2011-08-04-2/data007');
d08 = load_data('2011-08-04-2/data008');
d09 = load_data('2011-08-04-2/data009');
d10 = load_data('2011-08-04-2/data010');
d11 = load_data('2011-08-04-2/data011');

dataruns = {'d01str' 'd05str' 'd06' 'd07' 'd08' 'd09' 'd10' 'd11'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


% EIs
dataruns = {'d01str' 'd05str' 'd06' 'd07' 'd09' 'd10'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_ei(datarun, 'all');
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


% STAs and cones
dataruns = {'d01str' 'd05str' 'd08' 'd11'};
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = set_polarities(datarun);
    datarun = get_sta_fits_from_vision(datarun);
    
    datarun = load_cones(datarun, 1);
    datarun = make_mosaic_struct(datarun);
    datarun = make_voronoi_masks(datarun);
    
    assignin('caller', dataruns{i}, datarun);
end; clear dataruns datarun i


%% Picks
off_parasols_d01str = [14 27 32 38 44 50 97 165 201 234 289 244 297 266 313]; % we used these but the map moved
off_parasols_d05str = [317 2357 2836 4111 4279 4461 5779 6528 7216]; % Moved again!  A few may be usable.  Fuckers


%% Load stimuli
% Additity
stimuli = {'s06', 's10'};
mapdirs.s06 = '2011-08-04-2_f01_1and2';
mapdirs.s10 = '2011-08-04-2_f05_1and2';
for i = 1:length(stimuli)
    stim = read_stim_lisp_output(['2011-08-04-2/' stimuli{i}]);

    stim.mapdir = [server_data_path '2011-08-04-2/' mapdirs.(stimuli{i})];
    stim.mapims = {};
    for j = stim.mapindices
        fprintf('.');
        stim.mapims{j+1} = dlmread([stim.mapdir sprintf('/map-%.4d.txt', j)]);
    end; clear j mapims
    fprintf('\n');
    
    assignin('caller', stimuli{i}, stim);
end; clear stim stimuli i mapdirs

% Wack
stimuli = {'s07' 's09'};
mapdirs.s07 = '2011-08-04-2_f01_1and2all';
mapdirs.s09 = '2011-08-04-2_f05_1and2all';
for i = 1:length(stimuli)
    stim.mapdir = [server_data_path '2011-08-04-2/' mapdirs.(stimuli{i})];
    stim.mapims{1} = dlmread([stim.mapdir sprintf('/map-%.4d.txt', 0)]);
    stim.mapnyc{1} = map2manhattan(stim.mapims{1}, 'quiet', true);
    fprintf('.\n');
    assignin('caller', stimuli{i}, stim);
end; clear stim stimuli i mapdirs


%% Map d06 to d01str on EIs, get handpicked off parasols
cell_list_map = map_ei(d01str, d06, 'master_cell_type', {2});
off_parasols_fd01str_d06 = cell_list_map(get_cell_indices(d01str, off_parasols_d01str));
clear cell_list_map;

% The upshot
off_parasols_fd01str_d06 = {3860 4278 4340 [] 4728 4923 473 2163 5446 6891 2693 7205 2840 2358 []};


%% Map d07 to d01str on EIs, get handpicked off parasols
cell_list_map = map_ei(d01str, d07, 'master_cell_type', {2});
off_parasols_fd01str_d07 = cell_list_map(get_cell_indices(d01str, off_parasols_d01str));
clear cell_list_map;

% The upshot
off_parasols_fd01str_d07 = {3826 4282 4263 4323 4728 5056 918 2266 5446 6679 2958 7204 [] 2327 []};


%% Map d09 to d05str on EIs, get handpicked off parasols
cell_list_map = map_ei(d05str, d09, 'master_cell_type', {2});
off_parasols_fd05str_d09 = cell_list_map(get_cell_indices(d05str, off_parasols_d05str));
clear cell_list_map;

% The upshot
off_parasols_fd05str_d09 = {320 1593 [] 4043 4009 4325 5761 6411 7216};


%% Map d10 to d05str on EIs, get handpicked off parasols
cell_list_map = map_ei(d05str, d10, 'master_cell_type', {2});
off_parasols_fd05str_d10 = cell_list_map(get_cell_indices(d05str, off_parasols_d05str));
clear cell_list_map;

% The upshot
off_parasols_fd05str_d10 = {7368 1597 2823 4117 4233 4611 5868 6542 7351};


%% Build additivity struct
additivity_d06.conerun = d01str;
additivity_d06.conerun.addrgcs = off_parasols_d01str;
additivity_d06.addrun  = d06;
additivity_d06.addstim = s06;
additivity_d06.addrun.addrgcs  = off_parasols_fd01str_d06;
additivity_d06.stablerun = d08;


%% Additivity raster plot d06; need to check stability too though
additivity_compound_plot(additivity_d06, 'all');


%% Contrast response analysis d06
[cr06, templatex, templateyon, templateyoff] = additivity_contrast_response({d06 d10}, {s06 s10}, {off_parasols_fd01str_d06 off_parasols_fd05str_d10}, 0.75, 0.02, 'start', -0.1);
cr06(:,4,:) = cr06(:,1,:) + cr06(:,2,:);


%% Scatterplot d06
figure; hold on;
plot(squeeze(cr06(:,4,:)), squeeze(cr06(:,3,:)), 'b.');
plot(-10:350, -10:350, 'k--');
axis equal tight;
title 'Off Parasols 2011-08-04-2 d06';
xlabel 'Linear Prediction';
ylabel 'Observed Firing Rate';


%% Build additivity struct
additivity_d10.conerun = d05str;
additivity_d10.conerun.addrgcs = off_parasols_d05str;
additivity_d10.addrun  = d10;
additivity_d10.addstim = s10;
additivity_d10.addrun.addrgcs  = off_parasols_fd05str_d10;
additivity_d10.stablerun = d11;


%% Additivity raster plot d10; need to check stability too though
additivity_compound_plot(additivity_d10, 'all');
% Hard to tell, but I think some of these actually look okay


%% Contrast response analysis d10
[cr10, templatex, templateyon, templateyoff] = additivity_contrast_response({d10 d06}, {s10 s06}, {off_parasols_fd05str_d10 off_parasols_fd01str_d06}, 0.75, 0.02, 'start', -0.1);
cr10(:,4,:) = cr10(:,1,:) + cr10(:,2,:);


%% Scatterplot d10
figure; hold on;
plot(squeeze(cr10(:,4,:)), squeeze(cr10(:,3,:)), 'b.');
plot(-10:350, -10:350, 'k--');
axis equal tight;
title 'Off Parasols 2011-08-04-2 d10';
xlabel 'Linear Prediction';
ylabel 'Observed Firing Rate';


%% The wack, d07
s07.seeds = repmat(0:149, 1, 8);
s07 = get_wack_intensities(s07);

wack07.wackrun = d07;
wack07.wackrun.wackrgcs = off_parasols_fd01str_d07;
wack07.wackstim = s07;
wack07.conerun = d01str;
wack07.conerun.conergcs = off_parasols_d01str;
wack07.stablerun = d08;

figure; wack_compound_plot(wack07, 1:4);
figure; wack_compound_plot(wack07, 5:8);
figure; wack_compound_plot(wack07, 9:12);
figure; wack_compound_plot(wack07, 14);


% Both this and below are suspicious in that the cell with good
% cancellation has significantly higher baseline firing rate; suggests that
% it may not be an issue of same subunit so much as an issue of
% interneurons with high baselines not rectifying well.  Not clear why
% a particular RGC would have multiple interneurons with higher baseline
% though; it has to be both interneurons with high baseline to explain
% cancellation failure for both directions of up/down down/up.
%
% Nice example, except stability on 2nd sucks
%figure; wack_compound_plot(wack07, [3 8]);


%% The wack, d09
s09.seeds = repmat(0:149, 1, 8);
s09 = get_wack_intensities(s09);

wack09.wackrun = d09;
wack09.wackrun.wackrgcs = off_parasols_fd05str_d09;
wack09.wackstim = s09;
wack09.conerun = d05str;
wack09.conerun.conergcs = off_parasols_d05str;
wack09.stablerun = d11;

figure; wack_compound_plot(wack09, [1 2 4 5]); % 3 was not found
figure; wack_compound_plot(wack09, 6:9);

% See above for baseline firing issue
%
% Pretty good example
%figure; wack_compound_plot(wack09, [4 5 6 7]);