% Nishal P. Shah , August 2014

% Have a dataset ready with cells you want to nullify grouped under a
% common name(use Vision software).

startup_bertha

startup_rooster

datafile='nishal/2014-08-20-2/data001/data001';
type_name_inp='Good class';
no_images_per_movie=10;
start_image=10;
destination_mat='/Volumes/Analysis/nishal';
dest_raw='/Volumes/Data/stimuli/movies/null_space/';

solver=1; % 1 for LSQR  2 for CRAIG 
%%
[mov_orignial,mov_modify_new]=null_space_movie(datafile,type_name_inp,no_images_per_movie,start_image,destination_mat,solver);
% Doubt - do it have to 127.5 or not???

%% Original, stix 10
destination_raw_orignial=[dest_raw,'original_s10.rawMovie'];
write_movie([destination_mat,'/original/movie.mat'],destination_raw_orignial,10);
%% Modified, stix 10

destination_raw_modified=[dest_raw,'modified_s10.rawMovie'];
write_movie([destination_mat,'/modified/movie.mat'],destination_raw_modified,10);


%%
% write_movie('/Volumes/Analysis/nishal/sta_mov.mat','/Volumes/Data/stimuli/movies/null_space/sta.rawMovie');

%% Original, stix 2
movw=mov_orignial;
p3=zeros(size(movw,1)*5,size(movw,2)*5,size(movw,3));

for itime=1:size(movw,3)
p=movw(:,:,itime);
p2 = imresize(p,5,'nearest');
p3(:,:,itime)=double(imapprox(p2+127.5,gray,'dither'));
end

mov=p3 + 127.5;

save([destination_mat,'/original/mov_s2.mat'],'mov');

destination_raw_modified=[dest_raw,'original_d_s2.rawMovie'];
write_movie([destination_mat,'/original/mov_s2.mat'],destination_raw_modified,2);


%% Modified, stix 2
movw=mov_modify_new;
p3=zeros(size(movw,1)*5,size(movw,2)*5,size(movw,3));

for itime=1:size(movw,3)
p=movw(:,:,itime);
p2 = imresize(p,5,'nearest');
p3(:,:,itime)=double(imapprox(p2+127.5,gray,'dither'));
end

mov=p3 + 127.5;

save([destination_mat,'/modified/mov_s2.mat'],'mov');

destination_raw_modified=[dest_raw,'modified_d_s2.rawMovie'];
write_movie([destination_mat,'/modified/mov_s2.mat'],destination_raw_modified,2);


%% Modified, constrast enhance ?

mm_new = mov_modify_new*1.3;


mm_new(mm_new>127.5)=127.5;
mm_new(mm_new<-127.5)=-127.5;
mm_new=round(mm_new);

mov = mm_new+255/2;


save([destination_mat,'/modified/mov_con_s10.mat'],'mov');

destination_raw_modified=[dest_raw,'modified_con_s10.rawMovie'];
write_movie([destination_mat,'/modified/mov_con_s10.mat'],destination_raw_modified,10);
% last entry is stixel size

%% 
movie_white_len=120*100;%60*30;
earliest_peak_time=5; % select the most recent time in past where the STA value is high .. This gives us some stability in stimulus formation!!
[mov_modify_new]=white_noise(datafile,type_name_inp,destination_mat,movie_white_len,earliest_peak_time);