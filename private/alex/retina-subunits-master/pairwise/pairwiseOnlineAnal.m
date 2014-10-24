function pairwiseOnlineAnal(dat,CELLID)
%%
dat.cellNum = CELLID; %4621 and 2822 for 12-13-2011
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
fprintf('%-----------------------------%')
fprintf('Analyzing %s, #%g\r',dat.cellType{1},dat.rgcId);
fprintf('%-----------------------------%')

%%
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
 saveVmat(dat,fit,v,'vMats');
% loadVmat(dat,'vMats');
%%
clf
subplot(1,2,1);
plotPairwiseCell(dat,v);
subplot(1,2,2);
imagesc(v); axis square; colorbar;
set(gca,'CLim',[-max(v(:)) max(v(:))]);
h = title(sprintf('%s, #%g',dat.cellType{1},dat.rgcId));
set(h,'FontSize',20);
set(gcf,'Position',[763 507 676 299]);
set(gcf,'Color',[0.3 0.3 0.3]);
mkdir(dat.analPath);
set(gcf,'InvertHardcopy','off')
print('-dpng',fullfile(dat.analPath,sprintf('%g-vmat',CELLID)));
%%
selectCones(dat,triu(v),fullfile(dat.analPath,sprintf('%g-cone-selection.txt',CELLID)));
%%
clf
fit = getSTC(train);
plotSTC(dat.locs_c,fit,6);
%%
fprintf('%-----------------------------%\r')
fprintf('STA predictive accuracy\r')
fprintf('%-----------------------------%\r')

fit.g.type = 'logexp1';
fit.g.p = [0 1];
fit.D = 1;
fit.d_s = [0.6 0.1*ones(1,fit.D-1)];
fit = fitSTC(train,fit,0);
out1 = evalFitSTC(test,fit)
%%
fprintf('%-----------------------------%\r')
fprintf('STC predictive accuracy, 6 filters\r')
fprintf('%-----------------------------%\r')

fit.g.type = 'logexp1';
fit.g.p = [0 1];
fit.D = 6;
fit.d_s = [0.6 0.1*ones(1,fit.D-1)];
fit = fitSTC(train,fit,0);
out2 = evalFitSTC(test,fit)
%%
axes('position',[0.4 0.7 0.2 0.2]);
text(0,0.5,sprintf('STA accuracy, %g',out1.r2));
text(0,0.25,sprintf('STC accuracy, %g',out2.r2));
axis off
print('-dpng',fullfile(dat.analPath,sprintf('%g-stc',CELLID)));
