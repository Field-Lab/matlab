% Nishal P. Shah , September 2014

% Have a dataset ready with cells you want to nullify grouped under a
% common name(use Vision software).

startup_bertha

startup_rooster

datafile='2013-10-10-0/data000/data000';%'nishal/2014-08-20-2/data001/data001';

% type_name_inp = 'userCellList' for a list of cells.

no_images_per_movie=10;
start_image=10;
tag='test_stim';
destination_mat=['/Volumes/Analysis/nishal/',tag];
if(~exist(destination_mat,'dir'))
mkdir(destination_mat);
end

dest_raw='/Volumes/Data/stimuli/movies/null_space/';

movies=cell(20,1);
% 'solver' .. 1 for LSQR, 2 for CRAIG, 3 for Fourier 

%% BW - 10 sec movie cell type
global_vars2
var64=64;

cell_params=struct();
cell_params.type_name_inp='On Parasol';
cell_params.cell_list=[]; % if type_name_inp = 'userCellList' 
cell_params.use_fits=2;
cell_params.STAlen=14;

mov_params=struct();
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10;
mov_params.mean=0*255;
mov_params.deviation=0.36*255;
mov_params.scaling_loss=0.005; % a number in [0,1], fraction of values that is changed by scaling.

solver=3;
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
movies{1}=mov_orignial;
movies{2}=mov_modify_new;
mov_idx=1;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);
mov_idx=2;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);

%% BW - 10 sec user List
cell_params=struct();
cell_params.type_name_inp='userCellList';%'Off Parasol'
cell_params.cell_list=[301,1816,3604,3933]; % if type_name_inp = 'userCellList' 
cell_params.use_fits=2;
cell_params.STAlen=14;
% TODO : Add region specific cell list

mov_params=struct();
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10;
mov_params.mean=0.25*255;
mov_params.deviation=0.24*255;
mov_params.scaling_loss=0.02; % a number in [0,1], fraction of values that is changed by scaling.

solver=3;
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
movies{3}=mov_orignial;
movies{4}=mov_modify_new;
mov_idx=3;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);
mov_idx=4;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);
%% NSEM - 10 sec cell type

cell_params=struct();
cell_params.type_name_inp='Off Parasol';%'Off Parasol'
cell_params.cell_list=[]; % if type_name_inp = 'userCellList' 
cell_params.use_fits=2;
cell_params.STAlen=14;

mov_params=struct();
mov_params.mov_type='nsem';
mov_params.start_image=50;
mov_params.no_images_per_movie=10;
mov_params.mean=0.25*255; % mean on [0-255] scale.
mov_params.deviation_plus=0.73*255;
mov_params.deviation_minus=0.23*255;
mov_params.scaling_loss=0.02; % a number in [0,1], fraction of values that is changed by scaling.

solver=3; 
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
movies{5}=mov_orignial;
movies{6}=mov_modify_new;
mov_idx=5;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);
mov_idx=6;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);

%% NSEM - 10 sec - user List
cell_params=struct();
cell_params.type_name_inp='userCellList';%'Off Parasol'
cell_params.cell_list=[1218,2342,1880,3633,5689,6363,6784]; % if type_name_inp = 'userCellList' 
cell_params.use_fits=2;
cell_params.STAlen=14;

% TODO : Add region specific cell list
mov_params=struct();
mov_params.mov_type='nsem';
mov_params.start_image=50;
mov_params.no_images_per_movie=10;
mov_params.mean=0.25*255; % mean on [0-255] scale.
mov_params.deviation_plus=0.73*255;
mov_params.deviation_minus=0.23*255;
mov_params.scaling_loss=0.02; % a number in [0,1], fraction of values that is changed by scaling.

solver=3; 
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
movies{7}=mov_orignial;
movies{8}=mov_modify_new;
mov_idx=7;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);
mov_idx=8;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);



%% Correct means etc left ? 
movie_list= [1,2,4,5,6,8];
movie_full=zeros(size(movies{1},1),size(movies{1},2),size(movies{1},3)*length(movie_list));
icnt=1;
for imov=movie_list
movie_full(:,:,(icnt-1)*size(movies{imov},3)+1:icnt*size(movies{imov},3))=movies{imov};
icnt=icnt+1;
end
mov_idx=18;
write_movie_idx(destination_mat,movie_full,mov_idx);
display(sprintf('Movie Length %d',size(movie_full,3)));

figure;
for itime=1:10:size(movie_full,3)
imagesc(movie_full(:,:,itime));
colormap gray
colorbar
caxis([0,255]);
title(sprintf('Movie time %f s',itime/120));
pause(0.01);
end

%% BW - 15 min user List
cell_params=struct();
cell_params.type_name_inp='userCellList';%'Off Parasol'
cell_params.cell_list=[4864,1681,4141,4218,4817,4951,6856,4353,4351,4218]; % if type_name_inp = 'userCellList' 
cell_params.use_fits=1;

% TODO : Add region specific cell list

mov_params=struct();
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*60*15;
mov_params.mean=0.25*255;
mov_params.deviation=0.24*256;
mov_params.scaling_loss=0.02; % a number in [0,1], fraction of values that is changed by scaling.

solver=3;
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
movies{20}=mov_orignial;
movies{13}=mov_modify_new;
mov_idx=20;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);
mov_idx=13;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);