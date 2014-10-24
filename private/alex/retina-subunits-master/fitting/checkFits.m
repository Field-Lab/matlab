function out = checkFits(dat,opts)

% checkFits(dat,cellType,fitType)
%
% check whether fits exist for a data structure,
% cell type, and fit type

nCells = getCellTypeNum(dat.cellTypes,opts.cellType);
ids = getTypeIds(dat.cellTypes,opts.cellType);

out = [];
for cellNum=1:nCells
  celldat = [];
  celldat.cellInd = getCellInds(dat.rgcIds,ids(cellNum));
  celldat.rgcId = ids(cellNum);
  celldat.percent = opts.percent;
  celldat.dataSet = dat.dataSet;
  celldat.dataPath = dat.dataPath;
  celldat.analPath = dat.analPath;

  switch opts.fitType
    case {'greedy','greedy-local'}
      res = loadFit(celldat,opts.fitType,[]);
    case 'LN'
      res = loadFit(celldat,'greedy',[]);
  end

  if ~isempty(res) 
    if res.out_LN.r2 > opts.r2thresh
      out = cat(1,out,cellNum);
    else
      fprintf('(checkFits) r2 less than %g for cell %g\n',opts.r2thresh,celldat.rgcId);
    end
  end
end