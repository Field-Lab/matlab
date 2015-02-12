% Try LSQR
mov=mov_log{1};

addpath('../lsqrSOL/');
addpath('../craigSOL/');
b= -Ax(stas,mov,movie_time,n_cell);
alpha=0.1;
beta=0.5;
momentum=0.7;

%%
damp=0;
atol=10^-9;
btol=10^-9;
conlim=0; % Doubt!
itnlim=1000;
show=1;

[ x, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var ]...
   = lsqrSOL( movie_time*n_cell, movie_time*n_cell, @dual_AAtx, -2*b(:), damp, atol, btol, conlim, itnlim, show );

nu=reshape(x,[movie_time,n_cell]);
mov_modify_new=(-0.5)*(Atx(stas,nu,filt_dim1,filt_dim2,filt_len,movie_time,n_cell))+mov;

%%
damp=0;
atol=10^-8;
btol=10^-8;
conlim=1.0e+12; % Doubt!
itnlim=1000;
show=1;
tic;
[ x, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var ]...
   = lsqrSOL( movie_time*n_cell, movie_time*filt_dim1*filt_dim2, @dual_AAtx2, b(:), damp, atol, btol, conlim, itnlim, show );
toc;

mov_modify_new=reshape(x,[filt_dim1,filt_dim2,movie_time])+mov;
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell));

%% Craig's Method

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
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell));


