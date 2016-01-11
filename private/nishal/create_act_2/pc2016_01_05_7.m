
% Nishal P. Shah , September 2014

% Have a dataset ready with cells you want to nullify grouped under a
% common name(use Vision software).

startup_bertha

startup_rooster

datafile='2016-01-05-7/streamed/data000/data000';%'nishal/2014-08-20-2/data001/data001';

% type_name_inp = 'userCellList' for a list of cells.

no_images_per_movie=10;
start_image=10;
tag='pc2016_01_05_7_data000';
destination_mat=['/Volumes/Lab/Users/bhaishahster/',tag];
if(~exist(destination_mat,'dir'))
mkdir(destination_mat);
end

dest_raw='/Volumes/Data/stimuli/movies/null_space/';

movies=cell(20,1);
% 'solver' .. 1 for LSQR, 2 for CRAIG, 3 for Fourier 


%% Iterated Spatial nulling (high contrast) - has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='ON Parasol';%'userCellList';
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

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

%% low contrast WN, spatial null
cell_params=struct();
cell_params.type_name_inp='ON Parasol';%'userCellList';
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
mov_params.mov_type='bw-precomputed';
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


movies{3}=mov_orignial;
movies{4}=mov_modify_new;
mov_idx=3;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=4;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);

%% low contrast WN, spatial null ,selected ON parasol cells
cell_params=struct();
cell_params.type_name_inp='userCellList';
cell_params.cell_list=[871,1099,2072,2435,3017,3318,4203,4336,4941];
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


movies{5}=mov_orignial;
movies{6}=mov_modify_new;
mov_idx=5;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=6;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);



%% Iterated - contrast controlled Spatial nulling  - short -0.36
% has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='ON Parasol';%'userCellList';
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.36*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=22; % 22 for short variance controlled nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{7}=mov_orignial;
movies{8}=mov_modify_new;
mov_idx=7;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=8;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% Iterated - contrast controlled Spatial nulling  - short -0.36, ON selected
% has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='userCellList';
cell_params.cell_list=[871,1099,2072,2435,3017,3318,4203,4336,4941];
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
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.36*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=22; % 22 for short variance controlled nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{9}=mov_orignial;
movies{10}=mov_modify_new;
mov_idx=9;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=10;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);



%% Iterated - contrast controlled Spatial nulling  - short -0.24
% has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='ON Parasol';%'userCellList';
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
mov_params.mov_type='bw-precomputed';
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

solver=22; % 22 for short variance controlled nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{11}=mov_orignial;
movies{12}=mov_modify_new;
mov_idx=11;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=12;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% Iterated - contrast controlled Spatial nulling  - short -0.24, ON selected
% has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='userCellList';
cell_params.cell_list=[871,1099,2072,2435,3017,3318,4203,4336,4941];
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

solver=22; % 22 for short variance controlled nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{13}=mov_orignial;
movies{14}=mov_modify_new;
mov_idx=13;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=14;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% Iterated - contrast controlled Spatial nulling  - short -0.3
% has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='ON Parasol';%'userCellList';
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.3*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=22; % 22 for short variance controlled nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{15}=mov_orignial;
movies{16}=mov_modify_new;
mov_idx=15;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=16;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% Iterated - contrast controlled Spatial nulling  - short -0.3, ON selected
% has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='userCellList';
cell_params.cell_list=[871,1099,2072,2435,3017,3318,4203,4336,4941];
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
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.3*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=22; % 22 for short variance controlled nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{17}=mov_orignial;
movies{18}=mov_modify_new;
mov_idx=17;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=18;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);



%% Iterated - contrast controlled Spatial nulling  - short -0.18
% has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='ON Parasol';%'userCellList';
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.18*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=22; % 22 for short variance controlled nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{19}=mov_orignial;
movies{20}=mov_modify_new;
mov_idx=19;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=20;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% Iterated - contrast controlled Spatial nulling  - short -0.18, ON selected
% has interval parameter - remember to carefully decide the scaling (not stretching) when doing Iterated methods.
cell_params=struct();
cell_params.type_name_inp='userCellList';
cell_params.cell_list=[871,1099,2072,2435,3017,3318,4203,4336,4941];
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
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.18*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=22; % 22 for short variance controlled nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{21}=mov_orignial;
movies{22}=mov_modify_new;
mov_idx=21;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=22;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% WN+ Null, spatial

cell_params=struct();
cell_params.type_name_inp='ON Parasol';%'userCellList';
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
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval =2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=15; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{23}=mov_orignial;
movies{24}=mov_modify_new;
mov_idx=23;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=24;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);


%% WN+ Null, spatial, ON selected

cell_params=struct();
cell_params.type_name_inp='userCellList';
cell_params.cell_list=[871,1099,2072,2435,3017,3318,4203,4336,4941];
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
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval =2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=15; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

movies{25}=mov_orignial;
movies{26}=mov_modify_new;
mov_idx=25;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=26;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);

