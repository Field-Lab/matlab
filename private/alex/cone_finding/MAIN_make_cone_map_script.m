clear

%% Input parameters

piece = '2015-10-29-1';
run = 'data001';
online = false;  % look in online folder (/Volumes/Acquisition/date)
streamed = true;  %  look for streamed data in offline folder (/Volumes/Analysis/date/streamed)

magic_number = 20
inclusion_radius = 1.5

%% load data
'/Volumes/Acquisition/Analysis/2015-10-29-1/data004/data004'
path2data = '/Volumes/Acquisition/Analysis/2015-10-29-1/data004/data004';
datarun = load_data(path2data);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = load_cones(datarun, {});
% datarun = load_cones_ath(datarun,magic_number);
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
min_neighbor_dist = 1.5;
max_self_dist     = 4;
spaced = space_voronoi_masks(datarun, min_neighbor_dist, max_self_dist);
cone_map = index_masks(spaced, num2cell(indexes));

figure
imagesc(cone_map)

dlmwrite('/Volumes/Analysis/2015-10-29-1/stimuli/maps/map_data004.txt', cone_map, 'delimiter', '\t', 'newline', 'pc')
a = dlmread('/Volumes/Analysis/2015-10-29-1/stimuli/maps/map_data004.txt')

figure;
imagesc(a)

%% Visualize % for use - not necessary
figure; imagesc(cone_map); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); axis equal tight
figure; imagesc(cone_map); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); colormap gray; axis equal
figure; imagesc(cone_map > 0); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); colormap gray; axis equal


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


%%%END OF CODE FOR CELLS PICKING%%%%


%% Picking RGCs
plot_cone_mosaic(datarun, 'cone_colors', [0 0 0]);
datarun = get_sta_fits_from_vision(datarun);            

% Overlay rf fits and IDs
plot_rf_summaries(datarun, {4}, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')
plot_rf_summaries(datarun, {2}, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
% manual selection by cell index
plot_rf_summaries(datarun, [721 726 2809 3005 3183 4222 4758 5090], 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')

% plot original sta with voronoi cone regions individually for each cell
for rgc = [2806 3317 4353 4831 5222 6395 6890]
    figure
    plot_voronoi_over_rf(datarun, rgc, 'highlight_ranked_cones', [1 2], 'highlight_rgba', [1 0 0; 0 1 0]);
end

% % Plot DLL
% load([single_cone_path conesB '/results']);
% load([single_cone_path conesB '/parameters']);
% imagesc(matrix_scaled_up((norm_image(bcf.dll)-0.5).^0.4,3),...
%     'xdata',[1 datarunB.stimulus.field_width]-(1-bcf_params.kernel_spacing),'ydata',[1 datarunB.stimulus.field_height]-(1-bcf_params.kernel_spacing))
% axis image


% %% Look at cones in RGCs
% plot_rf
% 
% 
% %% Read saved map.txt file
% loaded_map = dlmread('/jacob/snle/data/2011-05-11-6/add_1_2_data002_600x600/map-0000.txt');
% imagesc(loaded_map);
% max(loaded_map(:))
% 


%% All cones

% give cell numbers for each cell cluster (copy from vision?)
clear offM
offM{1}=[6512 6182 2059 5794]; % data001
offM{1}=[2026 5795 6302 6511]; % data002

offM{1}=[721 726 2809 3005 3183 4222 4758 5090] % data003 for 2013-10-10-5


% offM{2}=[5422 5568 5569 5581 5941];
% offM{3}=[31 121 7681 7698];

newIndices=[];
clusterLength=zeros(1,size(offM,2));
figure
plot(datarun.cones.centers(:,1),datarun.cones.centers(:,2),'b+')% vertically flipped  
hold on

for cellClusters=1:size(offM,2)
    cells = offM{cellClusters};
    
    [mosaic_weights, selection, extras] = select_cone_weights(datarun, cells,...
        'thresh', 0,'radius', [0 inclusion_radius], 'polarity', 1,...
        'contiguity', true,'scale', 3.0);
    
    selection(excludes,:) = 0;
%     connectivity = mosaic_weights .* selection; % keep weights continuous valued
    
    % Check clusters visually
    
    all_COI=[]; %accumulates indices of all cones in cells of interest
    clear indices
    for i = 1:length(cells)
        indices{i} = find(selection(:,i));
        all_COI=[all_COI; indices{i}];
    end
    all_COI=unique(all_COI);
    
    % get convex hull of the region of interest and find all cones centers
    % within this hull
    a=convhull(datarun.cones.centers(all_COI,1),datarun.cones.centers(all_COI,2));
    
    
    IN = inpolygon(datarun.cones.centers(:,1),datarun.cones.centers(:,2),...
        datarun.cones.centers(all_COI(a),1),datarun.cones.centers(all_COI(a),2));
    
    %plot results
    plot(datarun.cones.centers(all_COI(a),1),datarun.cones.centers(all_COI(a),2),'r-',...
        datarun.cones.centers(IN>0,1),datarun.cones.centers(IN>0,2),'g*',...
        datarun.cones.centers(all_COI,1),datarun.cones.centers(all_COI,2),'m*')    
%     plot(datarun.cones.centers(IN>0,1),datarun.cones.centers(IN>0,2),'g+',...
%         datarun.cones.centers(all_COI,1),datarun.cones.centers(all_COI,2),'m+')
    
    % update indices: single array of all cones within the hull
    newIndices(cellClusters,1:nnz(IN))=find(IN);
    clusterLength(cellClusters)=sum(IN);
    
end

plot_rf_summaries(datarun, [721 726 2809 3005 3183 4222 4758 5090], 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')


maxIndices=size(newIndices,2);

for i=1:size(offM,2)
    if clusterLength(i)<maxIndices        
        centerCone=mean(datarun.cones.centers(newIndices(i,1:clusterLength(i)),:));
        distances=(datarun.cones.centers(:,1)-centerCone(1)).^2+(datarun.cones.centers(:,2)-centerCone(2)).^2;
        distances(newIndices(newIndices>0))=max(distances);
        [~,tmp]=sort(distances);
        tmp=tmp(1:(maxIndices-clusterLength(i)));
        newIndices(i,clusterLength(i)+1:end)=tmp;
        plot(datarun.cones.centers(tmp,1),datarun.cones.centers(tmp,2),'k*')
    end
    [~,tmp]=sort(datarun.cones.centers(newIndices(i,:),1));
    newIndices(i,:)=newIndices(i,tmp);    
    
    a=convhull(datarun.cones.centers(newIndices(i,:),1),datarun.cones.centers(newIndices(i,:),2));
    plot(datarun.cones.centers(newIndices(i,a),1),datarun.cones.centers(newIndices(i,a),2),'k-')
    
end

indices=cell(1,maxIndices);

for i=1:maxIndices
    indices{i}=newIndices(:,i);
end

cm = index_masks(spaced, indices);
figure; imagesc(cm); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);

mapdirname = [piece '/allconesOffmidgets'];
mkdir('/Volumes/Acquisition/maps/', mapdirname);
dlmwrite(['/Volumes/Acquisition/maps/' mapdirname '/map-0002.txt'], cm, 'delimiter', '\t', 'newline', 'pc');
save(['/Volumes/Acquisition/maps/' mapdirname '/cone_indices'],'indices')
% end of ath' code




%%



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
[sorted_cone_weights, sorted_cone_indices] = sort(connectivity);
for i = 1:size(sorted_cone_indices, 1)
    rowinds    = sorted_cone_indices(i,:);
    rowweights = sorted_cone_weights(i,:);
    maps{i} = rowinds(rowweights > 0);
end

% % Get all the cones into a vector, shuffle
% [r,c] = find(connectivity);
% selected_cones = unique(r);
% shuffled_selected_cones = selected_cones(randperm(length(selected_cones)));
% 
% % Save maps of individual serial cones
% m = 0;
% for c = shuffled_selected_cones
%     cone_map_c = index_masks(spaced, {shuffled_selected_cones(c)});
%     dlmwrite(sprintf('/snle/acquisition/maps/2011-06-24-6_f08_allcones/map-%.4d.txt', m), cone_map_c, 'delimiter', '\t', 'newline', 'pc');
%     m = m+1;
% end

mapdirname = '2012-04-13-1_f06_allcones';
mkdir('/snle/acquisition/maps/', mapdirname);

m = 0;
for i = 1:length(maps)
%     disp(num2str(i));
    
    if isempty(maps{i})
        continue;
    end
    
    cone_map_i = index_masks(spaced, {maps{i}});
    dlmwrite(sprintf('/snle/acquisition/maps/%s/map-%.4d.txt', mapdirname, m), cone_map_i, 'delimiter', '\t', 'newline', 'pc');
    
    m = m+1;
end


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