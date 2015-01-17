 ihbasprs.ncols = 12;  
 ihbasprs.hpeaks = [.1 .2];  
 ihbasprs.b = .01;  

 dt=.000833;
 [iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt)
 plot(iht*1000,ihbasis)
 %%
 addpath(genpath('~/Nishal/matlab/private/nora/'));
 alpha=0.01;
 beta=.000833*120;
 b=1;
 nfilters=12

 dt=0.00833;
 [basis, phi,tvec] = construct_cosine_basis3(alpha,beta,b,nfilters,dt)