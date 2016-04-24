
% Nishal P. Shah , September 2014

% Have a dataset ready with cells you want to nullify grouped under a
% common name(use Vision software).

startup_bertha

startup_rooster

datafile='2016-04-21-8/streamed/data018/data018';%'nishal/2014-08-20-2/data001/data001';

% type_name_inp = 'userCellList' for a list of cells.

no_images_per_movie=10;
start_image=10;
tag='pc2016_04_21_8_data018';
destination_mat=['/Volumes/Lab/Users/bhaishahster/',tag];
if(~exist(destination_mat,'dir'))
mkdir(destination_mat);
end

dest_raw='/Volumes/Data/stimuli/movies/null_space/';

movies=cell(20,1);
% 'solver' .. 1 for LSQR, 2 for CRAIG, 3 for Fourier 

%% low contrast WN, spatial null
cell_params=struct();
cell_params.type_name_inp='OFF midget';%'userCellList';
cell_params.cell_list=[];
cell_params.STAlen=14;
cell_params.sta_spatial=sprintf('%s/stas_spatial.mat',destination_mat);
cell_params.use_fits=2; % 2, 0,0,2
cell_params.sta_spatial_method=4;%1,2 ,3,4
% cell_params.sta_spatial_method = 1 for just using 4th frame, 2 is for fitting spatial STA. 
% Use cell_params.use_fits=2 (clipping) if cell_params.sta_spatial_method = 1 and 
% use cell_params.use_fits=0 (no processing of STA) if
% cell_params.sta_spatial_method = 2;
% STA spatial null Method 3 = low rank, 4 = average waveform and use it ..  

mov_params=struct();
mov_params.mov_type='bw';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.24*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=8; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
a = load('/Volumes/Analysis/2016-04-21-8/cone_data/manual/map_data015_manual_complete_178.txt');
mov_vor = mov_orignial;
moviexx=zeros(600,800,600);
for ivor=1:size(mov_vor,1)-1
    ivor
[r,c] = find(a==ivor); 
moviexx(r,c,:) = repelem(mov_vor(ivor,1,:),numel(r),numel(c),1);
end
movies{5}=moviexx;

mov_vor = mov_modify_new;
moviexx=zeros(600,800,600);
for ivor=1:size(mov_vor,1)-1
    ivor
[r,c] = find(a==ivor); 
moviexx(r,c,:) = repelem(mov_vor(ivor,1,:),numel(r),numel(c),1);
end
movies{6}=moviexx;

mov_params.stixel=1;
mov_idx=5;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=6;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);

%% low contrast WN, spatial null
cell_params=struct();
cell_params.type_name_inp='OFF parasol';%'userCellList';
cell_params.cell_list=[];
cell_params.STAlen=14;
cell_params.sta_spatial=sprintf('%s/stas_spatial.mat',destination_mat);
cell_params.use_fits=2; % 2, 0,0,2
cell_params.sta_spatial_method=4;%1,2 ,3,4
% cell_params.sta_spatial_method = 1 for just using 4th frame, 2 is for fitting spatial STA. 
% Use cell_params.use_fits=2 (clipping) if cell_params.sta_spatial_method = 1 and 
% use cell_params.use_fits=0 (no processing of STA) if
% cell_params.sta_spatial_method = 2;
% STA spatial null Method 3 = low rank, 4 = average waveform and use it ..  

mov_params=struct();
mov_params.mov_type='bw';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255; 
mov_params.deviation=0.24*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=8; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
 a = load('/Volumes/Analysis/2016-04-21-8/cone_data/manual/map_data015_manual_complete_178.txt');

mov_vor = mov_orignial;
moviexx=zeros(600,800,600);
for ivor=1:size(mov_vor,1)-1
    ivor
[r,c] = find(a==ivor); 
moviexx(r,c,:) = repelem(mov_vor(ivor,1,:),numel(r),numel(c),1);
end
movies{7}=moviexx;

mov_vor = mov_modify_new;
moviexx=zeros(600,800,600);
for ivor=1:size(mov_vor,1)-1
    ivor
[r,c] = find(a==ivor); 
moviexx(r,c,:) = repelem(mov_vor(ivor,1,:),numel(r),numel(c),1);
end
movies{8}=moviexx;

mov_params.stixel=1;
mov_idx=7;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=8;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% low contrast WN, spatial null
cell_params=struct();
cell_params.type_name_inp='ON parasol';%'userCellList';
cell_params.cell_list=[];
cell_params.STAlen=14;
cell_params.sta_spatial=sprintf('%s/stas_spatial.mat',destination_mat);
cell_params.use_fits=2; % 2, 0,0,2
cell_params.sta_spatial_method=4;%1,2 ,3,4
% cell_params.sta_spatial_method = 1 for just using 4th frame, 2 is for fitting spatial STA. 
% Use cell_params.use_fits=2 (clipping) if cell_params.sta_spatial_method = 1 and 
% use cell_params.use_fits=0 (no processing of STA) if
% cell_params.sta_spatial_method = 2;
% STA spatial null Method 3 = low rank, 4 = average waveform and use it ..  

mov_params=struct();
mov_params.mov_type='bw-precomputed2';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255; 
mov_params.deviation=0.24*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=8; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
 a = load('/Volumes/Analysis/2016-04-21-8/cone_data/manual/map_data015_manual_complete_178.txt');

mov_vor = mov_orignial;
moviexx=zeros(600,800,600);
for ivor=1:size(mov_vor,1)-1
    ivor
[r,c] = find(a==ivor); 
moviexx(r,c,:) = repelem(mov_vor(ivor,1,:),numel(r),numel(c),1);
end
movies{7}=moviexx;

mov_vor = mov_modify_new;
moviexx=zeros(600,800,600);
for ivor=1:size(mov_vor,1)-1
    ivor
[r,c] = find(a==ivor); 
moviexx(r,c,:) = repelem(mov_vor(ivor,1,:),numel(r),numel(c),1);
end
movies{8}=moviexx;

mov_params.stixel=1;
mov_idx=7;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=8;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
