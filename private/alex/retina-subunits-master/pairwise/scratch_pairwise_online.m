%%
dat = loadData('laptop','2012-04-13-1/data002');
%%

dat.cellType = 'off midget';
dat.cellNum = 1726; %4621 and 2822 for 12-13-2011
dat.trainFrac = 0.8;
dat.subDivide = 10;
dat.thresh = 2;
dat.percent = 0.2;
dat.threshType = 'peak';
dat.zscore = 1;
dat.tempFilter = 1;
[train test dat] = loadCellData(dat,2,0);
%%
%ind1 = 1; ind2 = 12;
%out = getPairwiseInt(train,test,ind1,ind2,fit,0)
%plotPairwiseInt(train,ind1,ind2,'gaussian',out,fit);
%%
clear fit
fit.g.type = 'exp';
fit.g.p = 0; 
if dat.stas.polarities{dat.cellInd} == -1
    fit.f.type = 'rectLinearNeg';
    fit.f.p = 0;
else
    fit.f.type = 'rectLinear';
    fit.f.p = 0;
end
v = getPairwiseIntAll(train,test,fit);
%%
%saveVmat(dat,fit,v,'vMats');
v = loadVmat(dat,'vMats');
%%
clf
subplot(1,2,1);
plotPairwiseCell(dat,v>20);

%print('-dpng','tmp.png');
%%
selectCones(dat,triu(v),'tmp.txt');
%%
fitstc = getSTC(train);
plotSTC(dat.locs_c,fitstc,6);
%%
fitstc.g.type = 'logexp1';
fitstc.g.p = [0 1];
fitstc.D = 1;
fitstc.d_s = [0.6 0.1*ones(1,fitstc.D-1)];
fitstc = fitSTC(train,fitstc,0);
out = evalFitSTC(test,fitstc)