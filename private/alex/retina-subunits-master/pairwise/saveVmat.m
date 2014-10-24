function saveVmat(celldat,fit,mat,folder,fileAppend)

if ~isdir(strcat(celldat.analPath,'/',folder,'/'));
  mkdir(strcat(celldat.analPath,'/',folder,'/'));
end

v.mat = mat;
v.f = fit.f;
v.g = fit.g;
v.percent = celldat.percent;
v.threshType = celldat.threshType;
v.zscore = celldat.zscore;
v.tempFilter = celldat.tempFilter;

fileName = sprintf('vmat-%g',celldat.rgcId);

if exist('fileAppend','var')
  fileName = strcat(fileName,sprintf('-%s',fileAppend));
end

save(fullfile(strcat(celldat.analPath,'/',folder,'/'),fileName),'v');

