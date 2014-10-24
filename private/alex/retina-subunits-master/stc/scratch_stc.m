% do an stc analysis
%load('~/Dropbox/Projects/retinaSubunits/data/examples/2008-08-27-5/data003/528-off midget.mat');
fit = getSTC(train);
plotSTC(cellDat.locs_c,fit,6);
%%
fit.g.type = 'logexp1';
fit.g.p = [0 1];
fit.D = 3;
fit.d_s = [0.6 0.1*ones(1,fit.D-1)];
fit = fitSTC(train,fit,0);
evalFitSTC(test,fit)