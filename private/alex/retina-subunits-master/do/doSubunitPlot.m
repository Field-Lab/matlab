function doSubunitPlot(celldat,dat,train,test,mode,reduced)

res = loadFit(celldat,mode,[]);
fits.fit_SUB = res.fit_SUB;
fits.fit_LN = res.fit_LN;

res_SING = loadFit(celldat,'singletons',[]);
fits.fit_SING = res_SING.fit_SUB;


if exist('reduced','var') && reduced
  plotModel2d_alt(fits,celldat,dat,train,test,mode);
else
  plotModel2d(fits,celldat,dat,train,test,mode);
end

set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'PaperOrientation','landscape');
savePath = sprintf('subunit-%s-%s-%s',mode,num2str(celldat.percent),celldat.cellType);
saveName = sprintf('%g-subunit-%s',celldat.rgcId,mode);
saveName = strcat(saveName,'-',num2str(celldat.percent));

if exist('reduced','var') && reduced
  saveName = strcat(saveName,'-reduced');
end
  
if ~isdir(fullfile(celldat.figurePath,savePath))
    mkdir(fullfile(celldat.figurePath,savePath))
end

print('-dpdf',fullfile(celldat.figurePath,savePath,strcat(saveName,'.pdf')));
