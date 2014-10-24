function res = loadFit(celldat,mode,thresh)

switch mode
  case 'singletons'
    subFolder = 'subunit';
    saveName = sprintf('%g-subunit-%s',celldat.rgcId,mode);
    %if dat.percent ~= 0.2
    saveName = strcat(saveName,'-',num2str(celldat.percent));
    %end
  case {'cliques','cliques-rand'}
    subFolder = 'subunit';
    saveName = sprintf('%g-subunit-%s-%g',celldat.rgcId,mode,thresh);
    saveName = strcat(saveName,'-',num2str(celldat.percent));
  case 'cliques-with-singletons'
    subFolder = 'subunit';
    saveName = sprintf('%g-subunit-%s-%g',celldat.rgcId,mode,thresh);
  case 'clusters-thresh'
    subFolder = 'subunit';
    saveName = sprintf('%g-subunit-%s-%g-%g',celldat.rgcId,mode,thresh(1),thresh(2));
  case 'apclusters'
    subFolder = 'subunit';
    saveName = sprintf('%g-subunit-%s',celldat.rgcId,mode);
  case {'greedy','greedy-local'}
    subFolder = 'subunit';
    saveName = sprintf('%g-subunit-%s',celldat.rgcId,mode);
    %if dat.percent ~= 0.2
    saveName = strcat(saveName,'-',num2str(celldat.percent));
    %end
  case 'LN'
    subFolder = 'LN';
    saveName = sprintf('%g-LN',celldat.rgcId);
  case 'STC'
    subFolder = 'STC';
    saveName = sprintf('%g-STC',celldat.rgcId);
end

if isfield(celldat,'simopts')
  saveName = strcat(saveName,'-',sprintf('d%g',celldat.simopts.train_dur),'-',num2str(celldat.simopts.iter),'-',celldat.simopts.mode);
end

if isfile(fullfile(celldat.analPath,subFolder,strcat(saveName,'-fit.mat')));
  load(fullfile(celldat.analPath,subFolder,strcat(saveName,'-fit.mat')));
else
  res = [];
  sprintf('(loadFit) cannot find fit for cell %g\n',celldat.rgcId);
end