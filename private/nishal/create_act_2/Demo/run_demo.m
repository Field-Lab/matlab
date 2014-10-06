
%% Load data

load('null_data.mat');
global_vars
addpath('./lsqrSOL/');
addpath('./craigSOL/');
b= -Ax(stas,mov,movie_time,n_cell);
alpha=0.1;
beta=0.5;
momentum=0.7;

%%  LSQR
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

    
%% CRAIG
atol=1.00e-06;%10^-8;
btol=1.00e-06;%10^-8;
conlim=1.00e+300; 
itnlim=1000;
show=1;
tic;
[ x, y, istop, itn, rnorm, Anorm, Acond, xnorm ]...
   = craigSOL( movie_time*n_cell, movie_time*filt_dim1*filt_dim2, @dual_AAtx2, b(:), atol, btol, conlim, itnlim, show );
toc;

mov_modify_new=reshape(x,[filt_dim1,filt_dim2,movie_time])+mov;
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell))
mov_modify_new = mov_modify_new*127.5/max(abs(mov_modify_new(:)));
    
