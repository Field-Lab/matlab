%% d02str: BW 3-4-0.48-22222-200x200 60.35 
off_midgets_d02str = [11 38 76 90 95 101 116 124 130 146 167 197 216 222];
d02str = load_data('2011-07-05-4/streamed/data002/data002');
d02str = load_params(d02str);
d02str = load_neurons(d02str);
d02str = load_ei(d02str, 'all');


%% d05: Voronoi additivity f02 15frameon=~250ms 750ms trials
d05 = load_data('2011-07-05-4/data005');
d05 = load_params(d05);
d05 = load_neurons(d05);
d05 = load_ei(d05, 'all');


%% stimuli
s05 = read_stim_lisp_output('2011-07-05-4/s05');
s05.mapdir = [server_data_path '2011-07-05-4/2011-07-05-4_f02_1and2'];

mapstimuli = {s05};
for s = 1:length(mapstimuli)
    stim = mapstimuli{s};
    mapims = {};
    for i = stim.mapindices
        fprintf('.');
        mapims{i+1} = dlmread([stim.mapdir sprintf('/map-%.4d.txt', i)]);
    end; 
    mapstimuli{s}.mapims = mapims;
    fprintf('\n');
end; clear s i stim mapims;

s05 = mapstimuli{1};
clear mapstimuli;


%% Map d05 to d02str on EIs, get handpicked off parasols
cell_list_map = map_ei(d02str, d05, 'master_cell_type', {4});
off_midgets_fd02str_d05 = cell_list_map(get_cell_indices(d02str, off_midgets_d02str));
clear cell_list_map;

% The upshot
off_midgets_fd02str_d05{7} = int32(2101);


%% Additivity raster plot; need to check stability too though
additivity_raster_plot(d05, s05, off_midgets_fd02str_d05);


%% Contrast response analysis
[cr05, templatex, templateyon, templateyoff] = additivity_contrast_response({d05}, {s05}, {off_midgets_fd02str_d05}, 0.75, 0.02, 'start', -0.1);
cr05(:,4,:) = cr05(:,1,:) + cr05(:,2,:);


%% Scatterplot
figure; hold on;
plot(squeeze(cr05(:,4,:)), squeeze(cr05(:,3,:)), 'b.');
plot(-10:350, -10:350, 'k--');
axis equal tight;
title 'Off Midgets 2011-07-05-4';
xlabel 'Linear Prediction';
ylabel 'Observed Firing Rate';

