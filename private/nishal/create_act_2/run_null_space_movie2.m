% Nishal P. Shah , September 2014

% Have a dataset ready with cells you want to nullify grouped under a
% common name(use Vision software).

startup_bertha

startup_rooster

datafile='2013-10-10-0/data000';%'nishal/2014-08-20-2/data001/data001';

% type_name_inp = 'userCellList' for a list of cells.

no_images_per_movie=10;
start_image=10;
tag='test-tag';
destination_mat=['/Volumes/Analysis/nishal/',tag];
if(~exist(destination_mat,'dir'))
mkdir(destination_mat);
end

dest_raw='/Volumes/Data/stimuli/movies/null_space/';

movies=cell(15,1);
% 'solver' .. 1 for LSQR, 2 for CRAIG, 3 for Fourier 
%%

cell_params=struct();
cell_params.type_name_inp='userCellList';%'Off Parasol'
cell_params.cell_list=[32,301,333,2808]; % if type_name_inp = 'userCellList' 

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
movies{1}=mov_orignial;
movies{2}=mov_modify_new;

%%
cell_params=struct();
cell_params.type_name_inp='On Parasol';
cell_params.cell_list=[]; % if type_name_inp = 'userCellList' 

mov_params=struct();
mov_params.mov_type='bw';
mov_params.movie_time=120*10;
mov_params.mean=0.25*255;
mov_params.deviation=0.24*255;
mov_params.scaling_loss=0.02; % a number in [0,1], fraction of values that is changed by scaling.

solver=3;
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
movies{3}=mov_orignial;
movies{4}=mov_modify_new;
%%
cell_params=struct();
cell_params.type_name_inp='userCellList';
cell_params.cell_list=[32,301,333,2808]; % if type_name_inp = 'userCellList' 

mov_params=struct();
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10%60*15;
mov_params.mean=0.25*255;
mov_params.deviation=0.24*255;
mov_params.scaling_loss=0.02; % a number in [0,1], fraction of values that is changed by scaling.

solver=3;
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);
movies{4}=mov_orignial;
movies{5}=mov_modify_new;


%% Correct means etc left ? 
for mov_idx=1:6
write_movie_idx(destination_mat,movies{mov_idx},mov_idx);

end

%[stim,height,width,header_size] = get_raw_movie('/Volumes/Analysis/nishal/test-tag/12.rawMovie',1440,1);



%  write_movie_idx('/Volumes/Analysis/nishal/test-tag/',double((stas{2}(:,:,3)+0.5)*255>129)*255,20);
%  close all
%  imagesc(double((stas{2}(:,:,3)+0.5)*255>129)*255);
%  colorbar
%  colormap gray
%  axis image
%  
% caxis([0,255]);
% 
% xx=double((stas{2}(:,:,3)+0.5)*255>129)*255;
% 
% mov=repmat(xx,[1,1,1440]);
%  write_movie_idx('/Volumes/Analysis/nishal/test-tag/',mov,20);
% 
% 














%%

[mov_orignial_four,mov_modify_new_four]=null_space_movie2(datafile,cell_params,mov_params,3);
[mov_orignial_lsqr,mov_modify_new_lsqr]=null_space_movie2(datafile,cell_params,mov_params,1);
%%
mov_modify_new_four=mov_modify_new_four/norm(mov_modify_new_four(:));
mov_modify_new_lsqr=mov_modify_new_lsqr/norm(mov_modify_new_lsqr(:));

figure;
for itime=1:1200
subplot(3,2,1)
imagesc(mov_orignial_four(:,:,itime));
colormap gray
axis image
colorbar

subplot(3,2,2)
imagesc(mov_modify_new_four(:,:,itime));
colormap gray
axis image
colorbar

subplot(3,2,3)
imagesc(mov_orignial_lsqr(:,:,itime));
colormap gray
axis image
colorbar

subplot(3,2,4)
imagesc(mov_modify_new_lsqr(:,:,itime));
colormap gray
axis image
colorbar


subplot(3,2,5)
imagesc(mov_orignial_lsqr(:,:,itime)-mov_orignial_four(:,:,itime));
colormap gray
axis image
colorbar

subplot(3,2,6)
imagesc(mov_modify_new_lsqr(:,:,itime)-mov_modify_new_four(:,:,itime));
colormap gray
axis image
colorbar
pause(1/120);
end

%%
figure;
for itime=110:1440
    itime
subplot(2,1,1);
imagesc(mov_orignial_nsem(:,:,itime));
caxis([min(mov_orignial_nsem(:)),max(mov_orignial_nsem(:))]);
colormap gray
axis image
colorbar

subplot(2,1,2);
imagesc(mov_modify_new_nsem(:,:,itime));
caxis([min(mov_orignial_nsem(:)),max(mov_orignial_nsem(:))]);
colormap gray
axis image
colorbar

pause(1/120)
end