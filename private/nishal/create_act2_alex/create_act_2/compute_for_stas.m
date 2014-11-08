function [mov_orig,mov_modify_new] = compute_for_stas(sstas)

global_vars
stas=sstas;
filt_len=30;
n_cell=length(sstas);
%%
no_mov=10;%no_images_per_movie;%10;
mov_ch=1;%30/no_mov;
mov_log=cell(mov_ch,1);
% See movie
start_image=1;

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

mov_log{imov_chunk}=mov;
end

filt_dim1=var64;
filt_dim2=32;
mov_log{1}=mov_log{1}*127.5/max(abs(mov_log{1}(:))); %contrast correction



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
solver_touse=1;
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

%% Verify computation
see_movie
mov_orig=mov_log{1};

end
