datarun = load_data('/Volumes/Acquisition/Analysis/2016-03-17-2/data001/data001');
movie_descr = 'BW-3-6-0.48-11111-265x200-60.35.xml';
cell_types = {1,2,3,4,5}; 
wid = 265;
hei = 200;
% load('/Volumes/Analysis/2016-02-17-4/data001_denoised_sta', 'denoised_sta'); % if present

% additionbal parameters
datarun.names.nickname = '';
datarun.piece.rig = 'A'; % 'A' or 'B'
datarun.piece.optical_path_direction = 'below'; %'below' or 'above'
datarun.piece.display = 'oled1'; % 'oled' or 'crt' with 1 for A, 2 for B


%% load stuff
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);
cell_inds = get_cell_indices(datarun, cell_types);

% denoise
% collect STAs
all_sta = zeros(hei,wid,length(datarun.cell_ids));
for i = 1:length(datarun.cell_ids)
    tmp = datarun.stas.stas{i};
    t = zeros(1,6);
    for j=1:6
        sta=squeeze(tmp(:,:,:,j));
        t(j) = max(abs(sta(:)));
    end
    [~, frame] = max(t);
    sta=tmp(:,:,:,frame);
    all_sta(:,:,i)=double(squeeze(sta));
end
noise_sta = mean(all_sta(:,:,cell_inds),3);
for i=1:length(datarun.cell_ids)
    tmp = datarun.stas.stas{i}-repmat(noise_sta,1,1,1,6);
    datarun.stas.stas{i} = tmp;
end

load('/Volumes/Analysis/2016-03-17-2/cone_data/data001/bcf.mat', 'bcf', 'bcf_params');

datarun.cones.centers = [];
datarun.cones.types = [];
choose_magic_number(datarun,bcf,bcf_params);


%%
magic_number = 20;
path2save = ['/Volumes/Analysis/2016-03-17-3/cone_data/data003/denoised_bayes-msf_', int2str(magic_number),'/'];
mkdir(path2save)
keep_indices = (bcf.all_added_cones(:,7) + magic_number*bcf.all_added_cones(:,8)) > 0;
cone_centers = bcf.all_added_cones(keep_indices,[4 5]);
datarun.cones.centers = cone_centers;
save([path2save, 'cones'], 'cone_centers', 'bcf', 'bcf_params');

%%
datarun = make_mosaic_struct(datarun);
datarun = make_voronoi_masks(datarun);
masks = datarun.cones.mosaic.voronoi_masks;

min_neighbor_dist = 2;
max_self_dist = 4;
spaced = space_voronoi_masks(datarun, min_neighbor_dist, max_self_dist);
cone_map = index_masks(spaced, num2cell(1:length(masks)));

figure
imagesc(cone_map)

dlmwrite(['/Volumes/Analysis/2016-03-17-3/cone_data/data003/denoised_bayes-msf_', int2str(magic_number),'/map_data003.txt'], cone_map, 'delimiter', '\t', 'newline', 'pc')

a = dlmread(['/Volumes/Analysis/2016-03-17-3/cone_data/data003/denoised_bayes-msf_', int2str(magic_number),'/map_data003.txt']);
figure;
imagesc(a)


%%

tmp = load('/Volumes/Analysis/2016-02-17-4/cone_data/data001/denoised_bayes-msf_70/cones.mat', 'cone_centers')
tmp1 = load('/Volumes/Analysis/2016-02-17-4/cone_data/data001/denoised_bayes-msf_120/cones.mat', 'cone_centers')

