global myPicture cone_list h i
global new_cone_list
h=[]; i=[];
cone_list=[];
myPicture=imread('/Volumes/Data/Auxiliary/2014-04-10-2/spot/Image1_edited.jpg');
myPicture=fliplr(myPicture(:,:,1));

hf=figure;
set(hf,'Toolbar','figure')
colormap gray

sbh=subplot(1,1,1);
imagesc(myPicture)

% crop and find all electrodes positions, scale figure

hcrop=uicontrol('style','pushbutton','Units','normalized','position',[0.73 0.945 0.12 0.03],'string','map electrodes','fontsize',16,...
    'callback','map_electrodes');

hclick=uicontrol('style','pushbutton','Units','normalized','position',[0.6 0.945 0.05 0.03],'string','cones','fontsize',16,...
    'callback','click_cones');

hdelete=uicontrol('style','pushbutton','Units','normalized','position',[0.5 0.945 0.05 0.03],'string','delete','fontsize',16,...
    'callback','delete_cones');

hcoord=uicontrol('style','pushbutton','Units','normalized','position',[0.3 0.945 0.12 0.03],'string','recalc coord','fontsize',16,...
    'callback','recalc_coord(cone_list)');


%%


magic_number = '10.00';
inclusion_radius = 1.5;

piece = '2014-06-04-7';
run = 'data002';
path2data=['/Volumes/Analysis/' piece '/' run '/' run];

datarun = load_data(path2data);

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta',[]);

datarun.cones.centers=new_cone_list;

datarun = make_mosaic_struct(datarun);
datarun = make_voronoi_masks(datarun);
masks = datarun.cones.mosaic.voronoi_masks;

excludes = [];

indexes = 1:length(masks);
indexes = setdiff(indexes, excludes);
cone_map = index_masks(masks, num2cell(indexes));


min_neighbor_dist = 2;
max_self_dist     = 3.5;
spaced = space_voronoi_masks(datarun, min_neighbor_dist, max_self_dist);
cone_map = index_masks(spaced, num2cell(indexes));
%dlmwrite('/snle/acquisition/maps/2011-07-05-2_f01_voronoicones/map-0000.txt', cone_map, 'delimiter', '\t', 'newline', 'pc')


%% Visualize % for use - not necessary
figure; imagesc(cone_map); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); axis equal tight
figure; imagesc(cone_map); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); colormap gray; axis equal
figure; imagesc(cone_map > 0); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale); colormap gray; axis equal



indices=cell(1,size(new_cone_list,1));

for i=1:size(new_cone_list,1)
    indices{i}=i;
end

cm = index_masks(spaced, indices);
figure; imagesc(cm); hold on; voronoi(datarun.cones.centers(:,1) .* datarun.cones.mosaic.voronoi_masks_scale, datarun.cones.centers(:,2) .* datarun.cones.mosaic.voronoi_masks_scale);


dirName = detect(find_cone_data(datarun), @(cd) (regexp(cd, magic_number)));
dirName = [path2data(1:end-15) dirName,'/maps/'];

mkdir(dirName);
dlmwrite([dirName '/map-0000.txt'], cm, 'delimiter', '\t', 'newline', 'pc');
coneIndices=indices;
coneCenters=datarun.cones.centers;
voronoiMasksScale=datarun.cones.mosaic.voronoi_masks_scale;
save([dirName '/index_info'],'coneIndices','cm','coneCenters','voronoiMasksScale')

% check saved info
load('/Volumes/Analysis/2014-04-10-2/2014-04-10-2_data000_data000-bayes-msf_20.00-all_BW-2-8/maps/index_info.mat')
figure; imagesc(cm); hold on; voronoi(coneCenters(:,1) .*voronoiMasksScale, coneCenters(:,2) .* voronoiMasksScale);

% read saved map and plot it 
loaded_map = dlmread('/Volumes/Analysis/2014-04-10-2/2014-04-10-2_data000_data000-bayes-msf_20.00-all_BW-2-8/maps/map-0000.txt');
imagesc(loaded_map);
max(loaded_map(:)) % total number of cones to stimulate
