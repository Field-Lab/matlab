function doSTCplot(celldat,train,test)

stc_tmp = getSTC(train,1);
plotSTC(celldat,stc_tmp,6);

% put the fit quality metrics on the plot
axes('Position',[0.8 0.025 0.2 0.2]);
axis off;

res = loadFit(celldat,'LN');
out_LN = res.out_LN;
res = loadFit(celldat,'STC');
out_STC = res.out_STC;

text(0,0,sprintf('R2, LN = %.03g',out_LN.r2),'FontSize',14);
text(0,0.15,sprintf('R2, STC = %.03g',out_STC.r2),'FontSize',14);

% save the figure
set(gcf, 'PaperPositionMode', 'auto');

savePath = sprintf('stc-%s',celldat.cellType);

if ~isdir(fullfile(celldat.figurePath,savePath))
  mkdir(fullfile(celldat.figurePath,savePath));
end
print('-dpdf',fullfile(celldat.figurePath,savePath,sprintf('%g-stc.pdf',celldat.rgcId)));
