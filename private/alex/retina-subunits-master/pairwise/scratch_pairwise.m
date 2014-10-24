% do an stc analysis
clear all
path = '~/Dropbox/Projects/retinaSubunits/code3/examples';
dataSet = '/2008-08-27-5/data003/';
%dataSet = '/2008-08-27-5/data003/';
file = '1112-off-midget.mat';
%file = '4456-off parasol.mat';
load(fullfile(path,dataSet,file));
%%
dat = loadData('laptop','2012-01-27-1/data004');
%%
fit.f.type = 'rectLinearNeg';
fit.f.p = 0;
fit.g.type = 'exp';
fit.g.p = [0 1];
out = getPairwiseInt(train,test,6,8,fit,0)
%v = getPairwiseIntAll(train,test,fit);
%%
plotPairwiseInt(train,3,6,'gaussian',out,fit);
%%
fit.f.type = 'exp';
if strcmpi(dat.cellType{1},'off midget') || ...
        strcmpi(dat.cellType{1},'off parasol')
    fit.g.type = 'rectLinearNeg';
else
    fit.g.type = 'rectLinear';
end
%v = getPairwiseIntAll(train,test,fit,dat.crossVal);
%saveVmat(dat,fit,v);

%%
plotPairwiseCell(cellDat,cellDat.v);
set(gcf,'Position',[763 507 437 299]);
h = title(sprintf('%s, #%g',cellDat.cellType,c ellDat.rgcId));
set(h,'FontSize',20);