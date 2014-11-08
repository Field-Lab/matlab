function [mov_orig,mov_modify_new]=null_space_movie(datafile,type_name_inp,no_images_per_movie,start_image,destination,solver_touse)
% Written by Nishal P. Shah
% Around August 2014

%% Start parallel pool
if(0)
if ((exist('matlabpool')==2) && (matlabpool('size') == 0))
  try
    matlabpool open
  catch me
    warning('Failed to open parallel sessions using matlabpool:\n  %s\n',...
        me.message);
  end
end
end

if (exist('parpool')==2)
  try
    if (isempty(gcp('nocreate')))
      parpool
    end
  catch me
    warning('Failed to open parallel pool using parpool:\n  %s\n',...
        me.message);
  end
end
%%
global_vars
%datafile = '2013-10-10-0/data000';
type_name= cell(1,1);
type_name{1}=type_name_inp;

datarun=load_data(datafile)
datarun=load_sta(datarun)
datarun=load_params(datarun)

%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!

matlab_cell_ids=get_cell_indices(datarun,type_name);
stas=datarun.stas.stas(matlab_cell_ids);
n_cell=length(stas);


%% Load STAS

stas_new=cell(length(stas),1);
for icell=1:length(stas)
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    stas_new{icell}=st_temp;
end
stas=stas_new;
filt_len=size(stas{1},4);

%% Load movies

no_mov=no_images_per_movie;%10;
mov_ch=1;%30/no_mov;
mov_log=cell(mov_ch,1);
% See movie


movie_time=120*no_mov;%size(movie,3);

for imov_chunk=1:mov_ch
    imov_chunk
movie_log=zeros(320,160,movie_time);
icnt=0;
for imov=start_image:no_mov+start_image-1
    imov
    icnt=icnt+1;
movv=load(sprintf('/Volumes/Data/stimuli/movies/eye-movement/NS_brownian/nishal_test_chunks/movie_chunk_%d',imov));
movv=movv.movie;
movie_log(:,:,(icnt-1)*120+1:icnt*120)=movv;
end
movie2=movie_log;
clear movie_log
var64=64;
mov=zeros(var64,32,movie_time);
for itime=1:movie_time
mov(:,:,itime)=(imresize(movie2(:,:,itime),[var64,32],'bilinear','Antialiasing',true)-127.5); % Doubt!!
end

% % Show movie
% figure
% for itime=1:movie_time
%     subplot(2,1,1);
%     imagesc(movie2(:,:,itime)');
%     caxis([0,255]);
%     axis image
%     %title(sprintf('Time: %d',itime));
%     colorbar
%     
%     subplot(2,1,2);
%    imagesc(mov(:,:,itime))
%    caxis([-128,128]);
%    colormap gray
%    axis image
%    %title(sprintf('Time: %d',itime));
%    colorbar
%    
%    pause(1/120);
% end
mov_log{imov_chunk}=mov;
end

filt_dim1=var64;
filt_dim2=32;
mov_log{1}=mov_log{1}*127.5/max(abs(mov_log{1}(:))); %contrast correction

%% Save original movie

orig_movie_path=sprintf([destination,'/original/']);
if(~exist(orig_movie_path,'dir'))
mkdir(orig_movie_path);
end
mov=mov_log{1}+127.5;
save(sprintf([orig_movie_path,'movie.mat']),'mov');

%%
% Try LSQR
imov=1;
mov=mov_log{imov};

addpath('../lsqrSOL/');
addpath('../craigSOL/');
b= -Ax(stas,mov,movie_time,n_cell);
alpha=0.1;
beta=0.5;
momentum=0.7;

%%  LSQR - Solve , change it to CRAIG when that starts working!
if(solver_touse==1)
damp=0;
atol=10^-6;
btol=10^-6;
conlim=1.0e+300; % Doubt!
itnlim=1000;
show=1;

tic;
[ x, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var1 ]...
   = lsqrSOL( movie_time*n_cell, movie_time*filt_dim1*filt_dim2, @dual_AAtx2, b(:), damp, atol, btol, conlim, itnlim, show );
toc;

mov_modify_new=reshape(x,[filt_dim1,filt_dim2,movie_time])+mov;
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell))
mov_modify_new = mov_modify_new*127.5/max(abs(mov_modify_new(:)));
else
    
atol=1.00e-06;%10^-8;
btol=1.00e-06;%10^-8;
conlim=1.00e+300; % Doubt!..   1.0e+12
itnlim=1000;
show=1;
tic;
[ x, y, istop, itn, rnorm, Anorm, Acond, xnorm ]...
   = craigSOL( movie_time*n_cell, movie_time*filt_dim1*filt_dim2, @dual_AAtx2, b(:), atol, btol, conlim, itnlim, show );
toc;

mov_modify_new=reshape(x,[filt_dim1,filt_dim2,movie_time])+mov;
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell))
mov_modify_new = mov_modify_new*127.5/max(abs(mov_modify_new(:)));
    
end 

%% Save modified movie 
mov=mov_modify_new+127.5;

modified_movie_path=sprintf([destination,'/modified/']);
if(~exist(modified_movie_path,'dir'))
mkdir(modified_movie_path);
end
save(sprintf([modified_movie_path,'movie.mat']),'mov');

%% Verify computation
see_movie
mov_orig=mov_log{1};

save(['/Volumes/Analysis/',datafile,'_null_mov.mat']);


end