%% streamed/data002
on_midgets_d02str = [  22   36   45  149  147   30  155  203  168  170   84]; % In order!
d02str = load_data('2011-05-11-6/streamed/data002/data002');
d02str = load_params(d02str);
d02str = load_sta(d02str,'load_sta',[]);
d02str = set_polarities(d02str);
d02str = load_cones(d02str, 'acquisition');
d02str = make_mosaic_struct(d02str);

%% streamed/data003
d03s = load_data('2011-05-11-6/streamed/data003/data003');
d03s = load_params(d03s);
d03s = load_neurons(d03s);
d03s = load_sta(d03s, 'load_sta', []);
d03s = get_sta_fits_from_vision(d03s);

d03s = set_polarities(d03s);
d03s = load_cones(d03s, 1);
d03s = make_mosaic_struct(d03s);
d03s = make_voronoi_masks(d03s);


%% data001
d01 = load_data('2011-05-11-6/data001');
d01 = load_params(d01);
d01 = load_neurons(d01);
d01 = load_ei(d01, 'all');
d01 = load_sta(d01, 'load_sta', []);
d01 = get_sta_fits_from_vision(d01);
d01 = set_polarities(d01);

%% data002
on_midgets_d02    = [4652 5477 5867 6421 6511 4879 6902 3272 7442 7486 1246]; % In order!

d02 = load_data('2011-05-11-6/data002');
d02 = load_params(d02);
d02 = load_neurons(d02);
d02 = load_ei(d02, 'all');
d02 = load_sta(d02, 'load_sta', []);
d02 = get_sta_fits_from_vision(d02);

d02 = set_polarities(d02);
d02 = load_cones(d02, 1);
d02 = make_mosaic_struct(d02);
d02 = make_voronoi_masks(d02);

%% data003
d03 = load_data('2011-05-11-6/data003');
d03 = load_params(d03);
d03 = load_neurons(d03);
d03 = load_ei(d03, 'all');
d03 = load_sta(d03, 'load_sta', []);
d03 = get_sta_fits_from_vision(d03);
d03 = set_polarities(d03);

%% data007
d07 = load_data('2011-05-11-6/data007');
d07 = load_params(d07);
d07 = load_neurons(d07);
d07 = load_ei(d07, 'all');


%% data008
d08 = load_data('2011-05-11-6/data008');
d08 = load_params(d08);
d08 = load_neurons(d08);
d08 = load_ei(d08, 'all');


%% Stimuli
s07 = read_stim_lisp_output('2011-05-11-6/stimuli/s07');
s08 = read_stim_lisp_output('2011-05-11-6/stimuli/s08');
s07.mapdir = [server_data_path '/2011-05-11-6/add_1_2_data002_600x600'];
s08.mapdir = [server_data_path '/2011-05-11-6/serial_1_2_data002_manmidgets'];

mapstimuli = {s07 s08};
for s = 1:length(mapstimuli)
    stim = mapstimuli{s};
    mapims = {};
    mapnyc = {};
    for i = stim.mapindices
        fprintf('.');
        mapims{i+1} = dlmread([stim.mapdir sprintf('/map-%.4d.txt', i)]);
        mapnyc{i+1} = map2manhattan(mapims{i+1}, 'quiet', true);
    end; 
    mapstimuli{s}.mapims = mapims;
    mapstimuli{s}.mapnyc = mapnyc;
    fprintf('\n');
end; clear s i stim mapims mapnyc;

s07 = mapstimuli{1};
s08 = mapstimuli{2};
clear mapstimuli;


%% Map d07 to d02 on EIs
cell_list_map = map_ei(d02, d07, 'master_cell_type', {3});
on_midgets_fd02_d07 = cell_list_map(get_cell_indices(d02, on_midgets_d02));
clear cell_list_map;

% The upshot
on_midgets_fd02_d07 = {4652 5476 5866 6421 6511 4818 6721 3272 7458 7486 1397};


%% Plot rasters d07: additivity; need to check stability
additivity_raster_plot(d07, s07, on_midgets_fd02_d07);
good = [1 2 3 6];
okay = [4 9];
bad = [7 8 11];
terrible = [5 10];


%% Searching d07

% ON Parasols
d07.mapd02 = map_ei(d02, d07);
additivity_raster_plot(d07, s07, [d07.mapd02{get_cell_indices(d02, {1})}], 'figure_by_cell_id', true);
% Semi-decent: 4381 4489 5071
% Very small, mostly decrement responses:  5776



%% Map d08 to d02 on EIs
cell_list_map = map_ei(d02, d08, 'master_cell_type', {3});
on_midgets_fd02_d08 = cell_list_map(get_cell_indices(d02, on_midgets_d02));
clear cell_list_map;

% The upshot
on_midgets_fd02_d08 = {4651 5477 5867 6422 6511 4877 6902 3272 7458 7576 1396};


%% Plot rasters d08: crosstalk control
crosstalk_raster_plot(d08, s08, on_midgets_fd02_d08);


%% Contrast response analysis
clear cr07;
use = [good okay];
[cr07, templatex, templateyon, templateyoff] = additivity_contrast_response({d07 d08}, {s07 s08}, {on_midgets_fd02_d07(use), on_midgets_fd02_d08(use)}, 0.75, 0.02, 'start', -0.1);
cr07(:,4,:) = cr07(:,1,:) + cr07(:,2,:);


%% Scatterplot
figure; hold on;
plot(squeeze(cr07(:,4,:)), squeeze(cr07(:,3,:)), '.');
plot(-10:100, -10:100, 'k--');
axis equal tight;

% Looks pretty bad.  
% Should probably rerun analysis using EI mapping instead of projection mapping.
% d07 d08 are analyzed without mapping, but need to check EI content.




%% data009
on_midgets_d09    = {4652 5342 5866 6421 [] 4818 6901 [] 7456 [] []}; % In order!
d09 = load_data('2011-05-11-6/data009');
d09 = load_params(d09);
d09 = load_sta(d09, 'load_sta', []);
d09 = get_sta_fits_from_vision(d09);
d09 = set_polarities(d09);
d09 = load_cones(d09, 1);
d09 = make_mosaic_struct(d09);
d09.cones.mosaic.voronoi_masks_scale = 2;


%% Cone projective field; searching

% First lets take a look through all the on midgets as a control
cell_list_map = map_ei(d02, d07, 'master_cell_type', {3});
all_on_midgets_fd02_d07 = cell_list_map(get_cell_indices(d02, {3}));
clear cell_list_map;
additivity_raster_plot(d07, s07, all_on_midgets_fd02_d07, 'intensities', [-0.5 0.5], 'figure_by_cell_id', true);

% Looks pretty good; can pretty easily pick out almost all the ones that
% were intentionally stimulated.

% A few others seem very marginally stimulated: 
%        4624, 4801, 5853
%in d02: 3048, 4802, 5854

% Some others look like they may have been very slightly repressed!  Surround??
%         3766, 3812, 4218, 6037, 6468, 7039, 7696
% in d02: 3721, 3843, 4622, 5722, 6468, 7191, 7697



% Okay, lets looks at On Parasols
cell_list_map = map_ei(d02, d07, 'master_cell_type', {1});
all_on_parasols_fd02_d07 = cell_list_map(get_cell_indices(d02, {1}));
clear cell_list_map;
additivity_raster_plot(d07, s07, all_on_parasols_fd02_d07, 'intensities', [-0.5 0.5], 'figure_by_cell_id', true);

% Not too many; only strong ones are on the double up / double down cases
% Prety good one: 5071 (d02: 5071)
% Fairly strong looking reverse: 5581 (d02: 5581)
% Okay: 4381, 5776 (d02: 4441, 5761)
% Weak: 1820, 3108, 4489, 7338 (d02: 1531, 3108, 4501, 7443)


% Off Parasols
cell_list_map = map_ei(d02, d07, 'master_cell_type', {2});
all_off_parasols_fd02_d07 = cell_list_map(get_cell_indices(d02, {2}));
clear cell_list_map;
additivity_raster_plot(d07, s07, all_off_parasols_fd02_d07, 'intensities', [-0.5 0.5], 'figure_by_cell_id', true);
% Some good ones, including single cone!
%      3395, 4324, 4673, 5014, 5270, 6907, 7341
% d02: 3241, 4398, 4669, 4787, 5672, 6936, 7458


% Off Midgets
cell_list_map = map_ei(d02, d07, 'master_cell_type', {4});
all_off_midgets_fd02_d07 = cell_list_map(get_cell_indices(d02, {4}));
clear cell_list_map;
additivity_raster_plot(d07, s07, all_off_midgets_fd02_d07, 'intensities', [-0.5 0.5], 'figure_by_cell_id', true);
% Nothing!!


%% ON Parasol rasters
additivity_raster_plot(d07, s07, [4381 5071 5776], 'intensities', [-0.5 0.5], 'figure_by_cell_id', true); % Better
s07.spec.index_map = s07.mapdir;
s07 = load_stimulus_maps(s07, '2011-05-11-6');

for id = [4441 5071 5761]
    figure
    plot_rf_stimmap(d02, id, s07.mapnycpoly{1}, 'colors', 'r', 'fit', false, 'color_transform', ones(3));
    bnds = axis;
    patchpolylines(s07.mapnycpoly{2}, 'colors', 'g', 'scale', [0.5 0.5])
    axis(bnds);
end


% SHOULD LOOK AT data008, the serial stimulation at high contrast.
s08.spec.index_map = s08.mapdir;
s08 = load_stimulus_maps(s08, '2011-05-11-6');
d08.mapd02 = map_ei(d02, d08);
d08.mapd07 = map_ei(d07, d08);

% 4381
crosstalk_raster_plot(d08, s08, d08.mapd07(get_cell_indices(d07, 4381)), 'maps', [0;11]);
figure(); plot_rf_stimmap(d02, 4441, s08.mapnycpoly{1}, 'autozoom', false, 'colors', 'y', 'fit', false, 'color_transform', ones(3));

% 5071
crosstalk_raster_plot(d08, s08, d08.mapd02(get_cell_indices(d02, 5071)));
figure(); plot_rf_stimmap(d02, 5071, s08.mapnycpoly{13}, 'autozoom', false, 'colors', 'y', 'fit', false, 'color_transform', ones(3));
% Hmm, doesn't look like much


% Try to check using the midgets
crosstalk_raster_plot(d08, s08, on_midgets_fd02_d08(2));
figure(); plot_rf_stimmap(d02, on_midgets_d02(2), s08.mapnycpoly{13}, 'autozoom', false, 'colors', 'y', 'fit', false, 'color_transform', ones(3));


% Combine STAs?
rf1 = get_rf(d01, 4441);
rf2 = get_rf(d02, 4441);
rf3 = get_rf(d03, 4308);
rfcompound = rf1+rf2+rf3;
datarun = d02;
datarun.stas.rfs{get_cell_indices(datarun, 4441)} = rfcompound;
figure
plot_rf_stimmap(datarun, 4441, s07.mapnycpoly{1}, 'colors', 'r', 'fit', false, 'color_transform', ones(3));
bnds = axis;
patchpolylines(s07.mapnycpoly{2}, 'colors', 'g', 'scale', [0.5 0.5])
axis(bnds);

rf1 = get_rf(d01, 5071);
rf2 = get_rf(d02, 5071);
rf3 = get_rf(d03, 5342);
rfcompound = rf1+rf2+rf3;
datarun = d02;
datarun.stas.rfs{get_cell_indices(datarun, 5071)} = rfcompound;
figure
plot_rf_stimmap(datarun, 5071, s07.mapnycpoly{1}, 'colors', 'r', 'fit', false, 'color_transform', ones(3));
bnds = axis;
patchpolylines(s07.mapnycpoly{2}, 'colors', 'g', 'scale', [0.5 0.5])
axis(bnds);

rf2 = get_rf(d02, 5761);
rf3 = get_rf(d03, 5761);
rfcompound = rf2+rf3;
datarun = d02;
datarun.stas.rfs{get_cell_indices(datarun, 5761)} = rfcompound;
figure
plot_rf_stimmap(datarun, 5761, s07.mapnycpoly{1}, 'colors', 'r', 'fit', false, 'color_transform', ones(3));
bnds = axis;
patchpolylines(s07.mapnycpoly{2}, 'colors', 'g', 'scale', [0.5 0.5])
axis(bnds);


%% Good RF On M
d01.mapd02 = map_ei(d02, d01);
d03.mapd02 = map_ei(d02, d03);

d02.onM = 5867;
d01.onM = d01.mapd02(get_cell_indices(d02, d02.onM));
d03.onM = d03.mapd02(get_cell_indices(d02, d02.onM));

d02.onMrf = get_rf(d02, d02.onM);
d01.onMrf = get_rf(d01, d01.onM{1});
d03.onMrf = get_rf(d03, d03.onM{1});
rfcompound = d01.onMrf + d02.onMrf + d03.onMrf;
datarun = d02;
datarun.stas.rfs{get_cell_indices(datarun, d02.onM)} = rfcompound;
plot_rf(datarun, d02.onM, 'color_transform', ones(3), 'fit', false, 'autozoom', true, 'scale', 10);
% All three together has nicest SNR.  Not quite as high a contrast as just d02+d03 but I think worth it.


%% Good RF On P
d01.mapd02 = map_ei(d02, d01);
d03.mapd02 = map_ei(d02, d03);

% 5761, 5929, 6361
d02.onP = 5761;
d01.onP = d01.mapd02(get_cell_indices(d02, d02.onP));
d03.onP = d03.mapd02(get_cell_indices(d02, d02.onP));

d02.onPrf = get_rf(d02, d02.onP);
d01.onPrf = get_rf(d01, d01.onP{1});
d03.onPrf = get_rf(d03, d03.onP{1});
rfcompound = d01.onPrf + d02.onPrf + d03.onPrf;
datarun = d02;
datarun.stas.rfs{get_cell_indices(datarun, d02.onP)} = rfcompound;
plot_rf(datarun, d02.onP, 'color_transform', ones(3), 'fit', false, 'autozoom', true, 'scale', 10);


%% Cone projective field; plotting
plot_rf_summaries(d02, on_midgets_d02, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k')
plot_rf_summaries(d02, [5071 4441 5761], 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(d02, [3241, 4398, 4669, 4787, 5672, 6936, 7458], 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')

%additivity_raster_plot(d07, s07, [4324, 4381, 4652], 'intensities', [-0.5 0.5], 'figure_by_cell_id', true);
additivity_raster_plot(d07, s07, [5776, 5270, 5866], 'intensities', [-0.5 0.5], 'figure_by_cell_id', true); % Better
cone2 = get_cones_by_weight_rank(d02, 5867, [2]);
figure(); plot_voronoi_over_rf(d02, 5867, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0])
figure(); plot_voronoi_over_rf(d02, 5672, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0])
figure(); plot_voronoi_over_rf(d02, 5761, 'highlight_cones', cone2, 'highlight_rgba', [0 1 0])


% Not working at the moment...
% additivity_d07_coneprojective.conerun = d02;
% additivity_d07_coneprojective.conerun.addrgcs = [on_midgets_d02      5071 5581 4441 5761];
% additivity_d07_coneprojective.addrun  = d07;
% additivity_d07_coneprojective.addstim = s07;
% additivity_d07_coneprojective.addrun.addrgcs  = [on_midgets_fd02_d07 5071 5581 4381 5776];
% additivity_d07_coneprojective.stablerun = d09;
% additivity_compound_plot(additivity_d07_coneprojective, 3);


%% Find specific matches from mapping
mapped = all_off_parasols_fd02_d07;
from = d02.cell_types{2}.cell_ids;
interest = [3395, 4324, 4673, 5014, 5270, 6907, 7341];

for i = 1:length(mapped)
    cell = mapped{i};
    if ~isempty(intersect(interest, cell))
        disp([num2str(cell) ': ' num2str(from(i))]);
    end
end

clear i cell mapped from interest ans

%% Old stuff below here
%% d07f02
d07f02 = load_data('2011-05-11-6/data007-from-data002');
d07f02 = load_neurons(d07f02);
d07f02 = load_ei(d07f02, handpicked_d02);


%% d08
d08 = load_data('2011-05-11-6/data008-from-data002');
d08 = load_neurons(d08);
d08 = load_ei(d08, handpicked_d02);


%% data009-from-data002
d09f02 = load_data('2011-05-11-6/data009-from-data002');
d09f02 = load_params(d09f02);
% d09f02 = load_sta(d09f02, 'load_sta', []);
% d09f02 = get_sta_fits_from_vision(d09f02);


%% Plot map with STA fits over
map = s08_mapim{9};
d = d02;
cells = handpicked_d02;

scale = 2;
ax1 = imagesc(map(1:scale:end, 1:scale:end)); colormap gray
hold on

plot_rf_summaries(d, cells, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b')
axis equal


%% Plot STA with voronoi over
d = d02;
cones = d02str.cones;
cells = handpicked_d02(9);
plot_sta(d, cells);
hold on;
voronoi(cones.centers(:,1), cones.centers(:,2));


%% Plot map with voronoi over
map = s08_mapim{9};
cones = d02str.cones;
voronoi_masks_scale = 2;
imagesc(map); hold on; voronoi(cones.centers(:,1) .* voronoi_masks_scale, cones.centers(:,2) .* voronoi_masks_scale);