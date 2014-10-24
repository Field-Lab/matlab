% do an stc analysis
clear all
path = '~/Dropbox/Projects/retinaSubunits/code3/examples';
dataSet = '/2010-03-05-2/data013/';
%dataSet = '/2008-08-27-5/data003/';
file = '677-on-midget.mat';
%file = '4456-off parasol.mat';
load(fullfile(path,dataSet,file));
%%
clear fit
fit.error = 'loglik';
fit.constraints.B_s = 'none';
fit.constraints.A_sc = 'L1';
fit.C = size(cellDat.locs_c,1);
fit.S = fit.C;
fit.locs_c = cellDat.locs_c;
fit.rf = cellDat.rf;
params.n = 7;
params.rand = 0;
params.thresh = 30;
params.res = 1;
init = 'clust';
switch init
    case {'singletons','doublets',...
            'doublets-overlap',...
            'singletons-and-doublets'};
        %fit.S = fit.C;
        init = mkSubAssgn(fit.C,init);
    case 'clust'
        init = getClusts('cliques',cellDat.v,params);
        init = init';
end
fit.I_sc = init;
fit.S = size(fit.I_sc,1);
%fit = addNeighbors(fit,0.4);
%%
fit.J_s = ones(1,size(fit.I_sc,1));
fit.f.type = 'rectLinear';
fit.f.p = [0 1];
fit.g.type = 'logexp1';
fit.g.p = [0 1];
fit.B_s = normMat(ones(1,fit.S),'L1',2);
fit.A_sc = normMat(rand(fit.S,fit.C).*fit.I_sc,'L1',2);
plotModel(fit);
%%
fit.f.type = 'spline';
fit.f.init = 'rectLinear';
fit.f.nknots = 8;
fit.f.smoothness = 3;
fit.f.extrap = 1;
fit.f.w = [];
fit = fitF(train,fit,1);
plotModel(fit);
%%
plotModel(fit); drawnow;
nIter =5;
for iIter = 1:nIter
    fit = fitB(train,fit,0);
    fit = fitA(train,fit,0);
    fit = fitF(train,fit,0);
    plotModel(fit); drawnow;
end