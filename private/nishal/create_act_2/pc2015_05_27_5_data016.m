
% Nishal P. Shah , September 2014

% Have a dataset ready with cells you want to nullify grouped under a
% common name(use Vision software).

startup_bertha

startup_rooster

datafile='2015-05-27-5/streamed/data016/data016';%'nishal/2014-08-20-2/data001/data001';

% type_name_inp = 'userCellList' for a list of cells.

no_images_per_movie=10;
start_image=10;
tag='pc2015_05_27_5_data016';
destination_mat=['/Volumes/Lab/Users/bhaishahster/',tag];
if(~exist(destination_mat,'dir'))
mkdir(destination_mat);
end

dest_raw='/Volumes/Data/stimuli/movies/null_space/';

movies=cell(20,1);
% 'solver' .. 1 for LSQR, 2 for CRAIG, 3 for Fourier 


%% Iterated Spatial nulling - has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='userCellList';%'userCellList';
cell_params.cell_list=[2191,3856,6166,1007,4926,1821,2311,5933,917,3137,5521,427,546,516,2268];
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/4;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-10-4-0.48-11111-32x32.xml';
mov_params.stixel=10;

mov_params.interval = 4; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=8; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{1}=mov_orignial;
movies{2}=mov_modify_new;
mov_idx=1;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=2;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);

%% WN+ Null, spatial

cell_params=struct();
cell_params.type_name_inp='userCellList';%'userCellList';
cell_params.cell_list=[2191,3856,6166,1007,4926,1821,2311,5933,917,3137,5521,427,546,516,2268];
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/4;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-10-4-0.48-11111-32x32.xml';
mov_params.stixel=10;

mov_params.interval =4; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=15; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{3}=mov_orignial;
movies{4}=mov_modify_new;
mov_idx=3;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=4;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);

%% Iterated Spatial nulling - has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='userCellList';%'userCellList';
cell_params.cell_list=[2191,3856,6166,1007,4926,1821,2311,5933,917,3137,5521,427,546,516,2268];
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/4;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-10-4-0.48-11111-32x32.xml';
mov_params.stixel=10;

mov_params.interval = 4; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=8; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_spacetag_movie2(datafile,cell_params,mov_params,solver);

movies{5}=mov_orignial;
movies{6}=mov_modify_new;
mov_idx=5;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=6;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);

%% WN+ Null, spatial

cell_params=struct();
cell_params.type_name_inp='userCellList';%'userCellList';
cell_params.cell_list=[2191,3856,6166,1007,4926,1821,2311,5933,917,3137,5521,427,546,516,2268];
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/4;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-10-4-0.48-11111-32x32.xml';
mov_params.stixel=10;

mov_params.interval =4; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=15; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{3}=mov_orignial;
movies{4}=mov_modify_new;
mov_idx=3;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=4;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% 
movie_list= [1,2,3,4,6,8];
movie_full=zeros(size(movies{1},1),size(movies{1},2),size(movies{1},3)*length(movie_list));
icnt=1;
for imov=movie_list
movie_full(:,:,(icnt-1)*size(movies{imov},3)+1:icnt*size(movies{imov},3))=movies{imov};
icnt=icnt+1;
end
mov_idx=19;
write_movie_idx(destination_mat,movie_full,mov_idx,mov_params.stixel);
display(sprintf('Movie Length %d',size(movie_full,3)));
%%
figure;
for itime=1:10:size(movie_full,3)
imagesc(movie_full(:,:,itime));
colormap gray
colorbar
caxis([0,255]);
title(sprintf('Movie time %f s',itime/120));
pause(0.01);
end


