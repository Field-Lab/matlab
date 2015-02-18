
% Nishal P. Shah , September 2014

% Have a dataset ready with cells you want to nullify grouped under a
% common name(use Vision software).

startup_bertha

startup_rooster

datafile='2011-09-02-0/data003/data003';%'nishal/2014-08-20-2/data001/data001';

% type_name_inp = 'userCellList' for a list of cells.

no_images_per_movie=10;
start_image=10;
tag='pc2014_11_05_9_data001_waste';
destination_mat=['/Volumes/Analysis/nishal/',tag];
if(~exist(destination_mat,'dir'))
mkdir(destination_mat);
end

dest_raw='/Volumes/Data/stimuli/movies/null_space/';

movies=cell(20,1);
% 'solver' .. 1 for LSQR, 2 for CRAIG, 3 for Fourier 



%% Spatio-temporal nulling. For different stixel sizes
cell_params=struct();
cell_params.type_name_inp='userCellList';%'OFF type 1';
cell_params.cell_list=[2314,1787,5597,5986,2087]; % if type_name_inp = 'userCellList' 
cell_params.use_fits=2;
cell_params.STAlen=14;

mov_params=struct();
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;
mov_params.scaling_loss=0.01; % a number in [0,1], fraction of values that is changed by scaling.
mov_params.stixel=4;
mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-4-1-0.48-11111-80x80.xml';


solver=5;
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
movies{1}=mov_orignial;
movies{2}=mov_modify_new;
mov_idx=1;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=2;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% Spatial nulling
cell_params=struct();
cell_params.type_name_inp='userCellList';%'OFF type 1';
cell_params.cell_list=[2314,1787,5597,5986,2087]; % if type_name_inp = 'userCellList'
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
mov_params.movie_time=120*10;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;
mov_params.scaling_loss=0.01; % a number in [0,1], fraction of values that is changed by scaling.
%mov_params.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-5-1-0.48-11111.xml'
mov_params.interval = 5; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.
mov_params.stixel=4;
mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-4-1-0.48-11111-80x80.xml';

solver=4; % Solver 4 used for spatial nulling!
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);


