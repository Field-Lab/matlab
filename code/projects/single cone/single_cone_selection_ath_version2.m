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
clear

magic_number = '140.00';
inclusion_radius = 1.5;

piece = '2014-09-10-1';
run = 'data001';
path2data=['/Volumes/Analysis/' piece '/Streamed/' run '/' run];

datarun = load_data(path2data);

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta',[]);
extra_mark=['all'];
datarun = load_cones(datarun,extra_mark);
% datarun = load_cones(datarun,magic_number);
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


%% Get cone weights for cone type of interest % for use

offM_d10=datarun.cell_types{1, 5}.cell_ids;


offM_d10=datarun.cell_types{1, 4}.cell_ids;

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, offM_d10,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', false,'scale', 3.0);
selection(excludes,:) = 0;
connectivity = mosaic_weights .* selection; % keep weights continuous valued
[~,sorted_cone_indices] = sort(connectivity);



%%%END OF CODE FOR CELLS PICKING%%%%


%% Picking RGCs
plot_cone_mosaic(datarun, 'cone_colors', [0 0 0]);
datarun = get_sta_fits_from_vision(datarun);            

% Overlay rf fits and IDs
plot_rf_summaries(datarun, {4}, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')
plot_rf_summaries(datarun, {2}, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
% manual selection by cell index
plot_rf_summaries(datarun, onM_d10, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')

offM{1}=offM_d10;

% plot original sta with voronoi cone regions individually for each cell
for rgc = offM{1}
    figure
    try
        plot_voronoi_over_rf(datarun, rgc, 'highlight_ranked_cones', [1 2], 'highlight_rgba', [1 0 0; 0 1 0]);
    catch err
    end
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
offM{1}=[346 812 826 1231 1606 1771 2131 2686 3032 3751 4006 5236 5342 5677 5882 6271 6621 7306];

offM{1}(3)=[];
offM{1}(offM{1}==2462)=[];
offM{1}(offM{1}==4036)=[];
offM{1}(offM{1}==4321)=[];


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
        'contiguity', false,'scale', 3.0);
    
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

plot_rf_summaries(datarun,offM{1}, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')


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



dirName = detect(find_cone_data(datarun), @(cd) (regexp(cd, magic_number)));
dirName = [path2data(1:end-15) dirName,'/maps/'];

mkdir(dirName);
dlmwrite([dirName '/map-0000.txt'], cm, 'delimiter', '\t', 'newline', 'pc');
coneIndices=indices;
cells=offM{:};
coneCenters=datarun.cones.centers;
voronoiMasksScale=datarun.cones.mosaic.voronoi_masks_scale;
save([dirName '/index_info'],'coneIndices','cells','cm','coneCenters','voronoiMasksScale')

% check saved info
load('/Volumes/Analysis/2014-04-10-2/2014-04-10-2_data000_data000-bayes-msf_20.00-all_BW-2-8/maps/index_info.mat')
figure; imagesc(cm); hold on; voronoi(coneCenters(:,1) .*voronoiMasksScale, coneCenters(:,2) .* voronoiMasksScale);

% read saved map and plot it 
loaded_map = dlmread('/Volumes/Analysis/2014-04-10-2/2014-04-10-2_data000_data000-bayes-msf_20.00-all_BW-2-8/maps/map-0000.txt');
imagesc(loaded_map);
max(loaded_map(:)) % total number of cones to stimulate
