% specify the directory to write to
% date_piece = '2010-10-18-1';
% datarun = 'data003';
% mosaic = 'on-midget';
% full_map_folder = 'all-cones';
% write_dir = '/snle/acquisition/maps/single-cone/';
% cd(write_dir)
% cone_folder = [date_piece,'-',datarun,'-',mosaic,];
% mkdir([date_piece,datarun,'-',full_map_folder])
% write_path = [write_dir, cone_folder,'/'];
% for use
piece = '2013-03-11-1';
run = 'data013';

% define data path
datarun = load_data(['/snle/acquisition/' piece '/' run '/' run]);
% path_and_name{1,2} = '_snle_acquisition_2011-01-19-0_data005_data005-bayes-msf_15.00-2011-01-19-0_data005';
% 
% datarun = load_data('/marte/snle/acquisition/2011-04-29-2/data004/data004');
% path_and_name{1,2} = '_snle_acquisition_2011-04-29-2_data004_data004-bayes-msf_20.00-2011-04-29-2_data004';
% datarun.cones.mosaic.voronoi_masks_scale = 2;

%datarun = load_data('2011-12-13-2/streamed/data018-0/data018-0');

% alternative data path
%datarun = load_data('2011-06-24-6/streamed/data003/data003');
%datarun.cones.mosaic.voronoi_masks_scale = 2;

% load data
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = load_cones(datarun); % load_cones(datarun, 'Analysis');
datarun = make_mosaic_struct(datarun);


%% make cone map new % for use
datarun = make_voronoi_masks(datarun);
masks = datarun.cones.mosaic.voronoi_masks;

excludes = [];

for i = 1:length(masks)
%     if sum(reshape(masks{i}, [], 1)) > 200
%         excludes(end+1) = i;
%     end
end; clear i

indexes = 1:length(masks);
indexes = setdiff(indexes, excludes);
cone_map = index_masks(masks, num2cell(indexes));


%% Spaced out % for use
min_neighbor_dist = 2;
max_self_dist     = 3.5;
spaced = space_voronoi_masks(datarun, min_neighbor_dist, max_self_dist);
cone_map = index_masks(spaced, num2cell(indexes));
%dlmwrite('/snle/acquisition/maps/2011-07-05-2_f01_voronoicones/map-0000.txt', cone_map, 'delimiter', '\t', 'newline', 'pc')


%% Visualize % for use - not necessary
figure; imagesc(cone_map); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); axis equal tight
figure; imagesc(cone_map); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); colormap gray; axis equal
figure; imagesc(cone_map > 0); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); colormap gray; axis equal


%% Voronoi cones %skip
maploc = '/snle/acquisition/maps/2011-12-13-2_f14_vorcones/';
dlmwrite([maploc 'map-0000.txt'], cone_map,  'delimiter', '\t', 'newline', 'pc')


%% Get cone weights for cone type of interest % for use
cells = [offM_d10];%array with cell ids selected in vision

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cells,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', false,'scale', 3.0);
selection(excludes,:) = 0;
connectivity = mosaic_weights .* selection; % keep weights continuous valued
[~,sorted_cone_indices] = sort(connectivity);
clear cells;


%% Strongest cone new %skip
peak_cone_indices = sorted_cone_indices(end,:);
cone_map1 = index_masks(spaced, {peak_cone_indices});


%% Second strongest new %skip
cone_map2 = index_masks(spaced, {sorted_cone_indices(end-1,:)});

%% Third strongest %skip
cone_map3 = index_masks(spaced, {sorted_cone_indices(end-2,:)});

%% Fourth strongest %skip
cone_map4 = index_masks(spaced, {sorted_cone_indices(end-3,:)});
maploc = '/snle/acquisition/maps/2011-12-04-1_f00_1234serial/';
dlmwrite([maploc 'map-0000.txt'], cone_map1,  'delimiter', '\t', 'newline', 'pc')
dlmwrite([maploc 'map-0001.txt'], cone_map2,  'delimiter', '\t', 'newline', 'pc')
dlmwrite([maploc 'map-0002.txt'], cone_map3,  'delimiter', '\t', 'newline', 'pc')
dlmwrite([maploc 'map-0003.txt'], cone_map4,  'delimiter', '\t', 'newline', 'pc')


%% 1st and 2nd new %skip
cone_map12 = index_masks(spaced, {sorted_cone_indices(end-1:end,:)});
maploc = '/snle/acquisition/maps/2011-08-04-6_f02_1and2/';
dlmwrite([maploc 'map-0000.txt'], cone_map1,  'delimiter', '\t', 'newline', 'pc')
dlmwrite([maploc 'map-0001.txt'], cone_map2,  'delimiter', '\t', 'newline', 'pc')
dlmwrite([maploc 'map-0002.txt'], cone_map12, 'delimiter', '\t', 'newline', 'pc')


%% Strongest four % for use. copy folder with results to the drive
cone_map1234 = index_masks(spaced, {sorted_cone_indices(end-3, :) sorted_cone_indices(end-2, :) sorted_cone_indices(end-1, :) sorted_cone_indices(end, :)});
figure; imagesc(cone_map1234); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); axis equal tight
maploc = ['/snle/acquisition/maps/' piece '/1234d10/'];
mkdir(maploc);
dlmwrite([maploc 'map-0000.txt'], cone_map1234,  'delimiter', '\t', 'newline', 'pc')




%%%END OF CODE FOR CELLS PICKING%%%%




%% Strongest four mixed with Freeman picks
picks = {sorted_cone_indices(end-3, :) sorted_cone_indices(end-2, :) sorted_cone_indices(end-1, :) sorted_cone_indices(end, :)};

% Initial for comparison
cone_map1234 = index_masks(spaced, picks);
figure; imagesc(cone_map1234); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); axis equal tight
lock

% Now substitute Jeremy's
for i = 1:length(jfcones_d05)
    cones = jfcones_d05{i};
    if isempty(cones), continue; end
    
    for j = 1:length(cones)
        picks{j}(i) = cones(j);
    end
end
clear i cones
cone_map1234 = index_masks(spaced, picks);
figure; imagesc(cone_map1234); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); axis equal tight


maploc = '/snle/acquisition/maps/2012-09-24-1/1234d03/';
dlmwrite([maploc 'map-0000.txt'], cone_map1234,  'delimiter', '\t', 'newline', 'pc');


%% Freeman four, new; cones just grouped into rows of 4 for each RGC
jfcones = jfcones_d13;

indexed_cones = mat2cell(jfcones', ones(size(jfcones,2),1))';
conemapjf4 = index_masks(spaced, indexed_cones);
figure; imagesc(conemapjf4); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);
maploc = ['/snle/acquisition/maps/' piece '/jffourd13/'];
mkdir(maploc);
dlmwrite([maploc 'map-0000.txt'], conemapjf4,  'delimiter', '\t', 'newline', 'pc')


%% Freeman four, old with cones already formatted for index_masks
% conemapjf4 = index_masks(spaced, offM_d04jf_cones);
% figure; imagesc(conemapjf4); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);
% maploc = '/snle/acquisition/maps/2012-09-06-0/jffourd04/';
% dlmwrite([maploc 'map-0000.txt'], conemapjf4,  'delimiter', '\t', 'newline', 'pc')


%% 1st and 2nd all different indices (for wacked 2; failed sparse attempt)
cone_map = index_masks(spaced, num2cell(sorted_cone_indices(end-1:end, :)));
dlmwrite('/snle/acquisition/maps/2011-08-04-2_f05_1and2all/map-0000.txt', cone_map, 'delimiter', '\t', 'newline', 'pc')


%% Strongest serial
%cone_map = index_masks(masks, num2cell(peak_cone_indices));


%% 1st and 2nd serial
%cone_map12serial = index_masks(spaced, num2cell(sorted_cone_indices(end-1:end,:)')); % Transpose for all 1st then all 2nd rather than cell 1 then cell 2 etc


%% 1st, 2nd, 3rd, 4th serial as series of maps
m = 0;
l = size(mosaic_weights, 1);
for i = [l-3 l-2 l-1 l]
    for j = 1:length(cells)
        cone_map_ij = index_masks(masks, {sorted_cone_indices(i, j)});
        dlmwrite(sprintf('/snle/acquisition/maps/2011-12-04-1_f00_1234serialsingle/map-%.4d.txt', m), cone_map_ij, 'delimiter', '\t', 'newline', 'pc');
        m = m+1;
    end
end


%% Strongest 2 in 2 serial groups
cone_map = index_masks(spaced, {sorted_cone_indices(end,:) sorted_cone_indices(end-1,:)});
dlmwrite('/snle/acquisition/maps2/wack/map-0000.txt', cone_map, 'delimiter', '\t', 'newline', 'pc')



%% Strongest 4 in 4 serial groups
cone_map = index_masks(spaced, {sorted_cone_indices(end,:) sorted_cone_indices(end-1,:) sorted_cone_indices(end-2,:) sorted_cone_indices(end-3,:)});
%write_cone_map = cone_map';
dlmwrite('/snle/acquisition/maps2/2011-01-19-0-data003.txt', write_cone_map, 'delimiter', '\t', 'newline', 'pc')


%% Strongest 4 in 4 separate files
cone_map1 = index_masks(spaced, {sorted_cone_indices(end,:)});
cone_map2 = index_masks(spaced, {sorted_cone_indices(end-1,:)});
cone_map3 = index_masks(spaced, {sorted_cone_indices(end-2,:)});
cone_map4 = index_masks(spaced, {sorted_cone_indices(end-3,:)});
dlmwrite('/snle/acquisition/maps2/map-0000.txt', cone_map1, 'delimiter', '\t', 'newline', 'pc')
dlmwrite('/snle/acquisition/maps2/map-0001.txt', cone_map2, 'delimiter', '\t', 'newline', 'pc')
dlmwrite('/snle/acquisition/maps2/map-0002.txt', cone_map3, 'delimiter', '\t', 'newline', 'pc')
dlmwrite('/snle/acquisition/maps2/map-0003.txt', cone_map4, 'delimiter', '\t', 'newline', 'pc')


%% Space out 2 strong cones per RGC (simulated annealing)
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cells,...
                                            'thresh', 0.5,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', false,'scale', 3.0);
selection(excludes,:) = 0;
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% starting point is just 1st and 2nd strongest as above
[~,sorted_cone_indices] = sort(connectivity);
initial = sorted_cone_indices(end-1:end,:)';

% Loss and generator functions for simulated annealing
loss = @(indexes) (cone_spacing_cost_function(datarun.cones.centers(indexes(:),:)));
generator = @(indexes) (rand_neighbor_cone_selection(indexes, selection));

% Get default annealing options, set generator and verbosity
saoptions = anneal();
saoptions.Generator = generator;
saoptions.Verbosity = 2; % Verbose

samin = anneal(loss, initial, saoptions);
cone_map_sa = index_masks(spaced, {samin(:)});


%% Picking RGCs
plot_cone_mosaic(datarun, 'cone_colors', [0 0 0]);
datarun = get_sta_fits_from_vision(datarun);            

% Overlay rf fits and IDs
plot_rf_summaries(datarun, {4}, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')
plot_rf_summaries(datarun, {2}, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarun, offM_d02, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarun, offM_d02, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')
plot_rf_summaries(datarun, [1298 393 1204 4156 7669 2133 4507 7506], 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarun, [181 619 1411 1966 1982 3422 4066 4489 1187 7760 6814 6093 4866 1096 5436 3258 2193], 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')


for rgc = [61 293 1237 1561 7007]
    figure
    plot_voronoi_over_rf(datarun, rgc, 'highlight_ranked_cones', [1 2], 'highlight_rgba', [1 0 0; 0 1 0]);
end

% Plot DLL
load([single_cone_path conesB '/results']);
load([single_cone_path conesB '/parameters']);
imagesc(matrix_scaled_up((norm_image(bcf.dll)-0.5).^0.4,3),...
    'xdata',[1 datarunB.stimulus.field_width]-(1-bcf_params.kernel_spacing),'ydata',[1 datarunB.stimulus.field_height]-(1-bcf_params.kernel_spacing))
axis image


%% Look at cones in RGCs
plot_rf


%% Read saved map.txt file
loaded_map = dlmread('/jacob/snle/data/2011-05-11-6/add_1_2_data002_600x600/map-0000.txt');
imagesc(loaded_map);
max(loaded_map(:))


%% Exclude cones manually?
datarun = get_sta_fits_from_vision(datarun);
for rgcid = onM_d04
    figure
    plot_rf(datarun, rgcid);
    hold on;
    plot_cone_mosaic(datarun, 'fig_or_axes', gca, 'clear', false, 'label', true);
    autozoom_to_fit(datarun, rgcid);
end


%% All cones

%Select center cones for RGCs of interest
cells = [off_midget_list_data003_additivity];

% extract connectivity
% [mosaic_weights, selection, extras] = select_cone_weights(datarun, cells,...
%     'thresh', 0,...
%     'radius', [0 2], 'polarity', 0,...
%     'contiguity', true, 'scale', 3.0);

[mosaic_weights, selection, extras] = select_cone_weights(datarun, cells,...
    'thresh', 0,...
    'radius', [0 1.5], 'polarity', 1,...
    'contiguity', true,'scale', 3.0);

selection(excludes,:) = 0;
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% Check clusters visually
for i = 1:length(cells)
    indices{i} = find(connectivity(:,i));
end
max(cellfun(@length, indices)) % How many?
cm = index_masks(spaced, indices);
figure; imagesc(cm); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);


[sorted_cone_weights, sorted_cone_indices] = sort(connectivity, 1, 'descend');
for i = 1:size(sorted_cone_indices, 1)
    rowinds    = sorted_cone_indices(i,:);
    rowweights = sorted_cone_weights(i,:);
    pixind = rowinds(rowweights > 0);
    
    if isempty(pixind), break; end
    pixindices{i} = pixind;
end


cm = index_masks(spaced, pixindices);
figure; imagesc(cm); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);

mapdirname = [piece '/allcones03'];
mkdir('/snle/acquisition/maps/', mapdirname);
dlmwrite(['/snle/acquisition/maps/' mapdirname '/map-0000.txt'], cm, 'delimiter', '\t', 'newline', 'pc');


%% All cones; midgets and parasols together!

M = offM_d09;
P = offP_d09;

% Midgets
[mosaic_weights, selection, extras] = select_cone_weights(datarun, M,...
    'thresh', 0.25,...
    'radius', [0 inf], 'polarity', 1,...
    'contiguity', true,'scale', 3.0);

selection(excludes,:) = 0;
Mconnectivity = mosaic_weights .* selection; % keep weights continuous valued

% Check clusters visually
for i = 1:length(M)
    Mindices{i} = find(Mconnectivity(:,i));
end
max(cellfun(@length, Mindices)) % How many?
cmM = index_masks(spaced, Mindices);
figure; imagesc(cmM); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);


% Parasols
[mosaic_weights, selection, extras] = select_cone_weights(datarun, P,...
    'thresh', 0.5,...
    'radius', [0 inf], 'polarity', 1,...
    'contiguity', true,'scale', 3.0);

selection(excludes,:) = 0;
Pconnectivity = mosaic_weights .* selection; % keep weights continuous valued

% Check clusters visually
for i = 1:length(P)
    Pindices{i} = find(Pconnectivity(:,i));
end
max(cellfun(@length, Pindices)) % How many?
cmP = index_masks(spaced, Pindices);
figure; imagesc(cmP); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);


[Psorted_cone_weights, Psorted_cone_indices] = sort(Pconnectivity, 1, 'descend');
for i = 1:size(Psorted_cone_indices, 1)
    rowinds    = Psorted_cone_indices(i,:);
    rowweights = Psorted_cone_weights(i,:);
    pixind = rowinds(rowweights > 0);
    
    if isempty(pixind), break; end
    Ppixindices{i} = pixind;
end


[Msorted_cone_weights, Msorted_cone_indices] = sort(Mconnectivity, 1, 'descend');
for i = 1:size(Msorted_cone_indices, 1)
    rowinds    = Msorted_cone_indices(i,:);
    rowweights = Msorted_cone_weights(i,:);
    pixind = rowinds(rowweights > 0);
    
    if isempty(pixind), break; end
    Mpixindices{i} = pixind;
end

% Put them together, reiterating the midgets during the parasols
% Doesn't work; a pixel can't have multiple map indices; would have to use
% multiple maps...
% MintoP = floor(length(Ppixindices) / length(Mpixindices));
% pixindices = Ppixindices;
% for i = 1:MintoP
%     for j = 1:length(Mpixindices);
%         ind = (i-1)*length(Mpixindices) + j;
%         pixindices{ind} = [pixindices{ind} Mpixindices{j}];
%     end
% end


% Put them together the boring way.  Add midgets in in reverse to balance
% things out though at least.
pixindices = Ppixindices;
for j = 1:length(Mpixindices);
    ind = length(pixindices) - j + 1;
    pixindices{ind} = [pixindices{ind} Mpixindices{j}];
end

cm = index_masks(spaced, pixindices);
figure; imagesc(cm); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);
mapdirname = '2012-09-13-2/allconesd09';
mkdir('/snle/acquisition/maps/', mapdirname);
dlmwrite(['/snle/acquisition/maps/' mapdirname '/map-0000.txt'], cm, 'delimiter', '\t', 'newline', 'pc');


%% Null stimulus all cones
spaced1 = space_voronoi_masks(datarun, min_neighbor_dist, 4);
cm1 = index_masks(spaced1, {indexes});

spaced2 = space_voronoi_masks(datarun, min_neighbor_dist, 6);
cm2 = index_masks(spaced2, {indexes});

cm3 = cm2.*2;
cm3(cm1 > 0) = 1;

sum(cm1(:))
sum(cm3(:) == 2)
figure; imagesc(cm3); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); colormap gray; axis equal

dlmwrite('/snle/acquisition/maps/2012-09-13-2/nulld05/map-0000.txt', cm3, 'delimiter', '\t', 'newline', 'pc')


%% Old multimap all cones (steps before "check clusters visually" are the same as above)
% [sorted_cone_weights, sorted_cone_indices] = sort(connectivity);
% for i = 1:size(sorted_cone_indices, 1)
%     rowinds    = sorted_cone_indices(i,:);
%     rowweights = sorted_cone_weights(i,:);
%     maps{i} = rowinds(rowweights > 0);
% end
% 
% % % Get all the cones into a vector, shuffle
% % [r,c] = find(connectivity);
% % selected_cones = unique(r);
% % shuffled_selected_cones = selected_cones(randperm(length(selected_cones)));
% % 
% % % Save maps of individual serial cones
% % m = 0;
% % for c = shuffled_selected_cones
% %     cone_map_c = index_masks(spaced, {shuffled_selected_cones(c)});
% %     dlmwrite(sprintf('/snle/acquisition/maps/2011-06-24-6_f08_allcones/map-%.4d.txt', m), cone_map_c, 'delimiter', '\t', 'newline', 'pc');
% %     m = m+1;
% % end
% 
% mapdirname = '2012-04-13-1_f06_allcones';
% mkdir('/snle/acquisition/maps/', mapdirname);
% 
% m = 0;
% for i = 1:length(maps)
% %     disp(num2str(i));
%     
%     if isempty(maps{i})
%         continue;
%     end
%     
%     cone_map_i = index_masks(spaced, {maps{i}});
%     dlmwrite(sprintf('/snle/acquisition/maps/%s/map-%.4d.txt', mapdirname, m), cone_map_i, 'delimiter', '\t', 'newline', 'pc');
%     
%     m = m+1;
% end


%% Null stimulus particular RGCs
name = 'nullf13';

%Select center cones for RGCs of interest
cells = [nullf13];

[mosaic_weights, selection, extras] = select_cone_weights(datarun, cells,...
    'thresh', 0.1,...
    'radius', [0 inf], 'polarity', 1,...
    'contiguity', true,'scale', 3.0);

selection(excludes,:) = 0;
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% Check clusters visually
for i = 1:length(cells)
    indices{i} = find(connectivity(:,i));
end
cm = index_masks(spaced, indices);

% Buffer around center
buffer = max_self_dist+0.75;
spaced2 = space_voronoi_masks(datarun, 0, buffer);
cm2 = index_masks(spaced2, indices);

% Filling area but not too close to other cones
ringwidth = 2.75;
spaced3 = space_voronoi_masks(datarun, buffer, buffer+ringwidth);
cm3 = index_masks(spaced3, indices);

% Buffered rings
cm4 = cm3;
cm4(cm2 > 0) = 0;
cm4(cm4 > 0) = 1;

% Centers
cm5 = cm;
cm5(cm5 > 0) = 1;

% Combine and inspect
cm6 = cm4.*2 + cm5;

% Diagnostics
figure; imagesc(cm6); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); colormap gray; axis equal tight
sum(cm4(:))
sum(cm5(:))

dirname = fullfile('/snle/acquisition/maps', piece, name);
mkdir(dirname);
dlmwrite(fullfile(dirname, 'map-0000.txt'), cm6, 'delimiter', '\t', 'newline', 'pc')