function [celldat,fit] = fitToJSON_load_cell(dat,opts,cellNum)

cellNum = opts.cellInds(cellNum);

ids = getTypeIds(dat.cellTypes,opts.cellType);

celldat = [];
celldat.cellInd = getCellInds(dat.rgcIds,ids(cellNum));
celldat.rgcId = ids(cellNum);
celldat.percent = opts.percent;
celldat.threshType = 'peak';
celldat.dataSet = dat.dataSet;
celldat.dataPath = dat.dataPath;
celldat.analPath = dat.analPath;
celldat.figurePath = dat.figurePath;
celldat.computer = dat.computer;
celldat.locations = dat.locations;

celldat = getConeInds(dat,celldat);

celldat.rf = dat.rfs{celldat.cellInd};
celldat.xRange = ([floor(min(celldat.locs_c(:,1))-1.5) ...
  ceil(max(celldat.locs_c(:,1))+1.5)]);
celldat.yRange = ([floor(min(celldat.locs_c(:,2))-1.5) ...
  ceil(max(celldat.locs_c(:,2))+1.5)]);

celldat.xRange(1) = max(celldat.xRange(1),ceil(min(dat.locations(:,1))));
celldat.xRange(2) = min(celldat.xRange(2),floor(max(dat.locations(:,1))));
celldat.yRange(1) = max(celldat.yRange(1),ceil(min(dat.locations(:,2))));
celldat.yRange(2) = min(celldat.yRange(2),floor(max(dat.locations(:,2))));

if ~isempty(celldat.rf)
  celldat.rf = celldat.rf(celldat.yRange(1):celldat.yRange(2),...
    celldat.xRange(1):celldat.xRange(2),:);
end

switch opts.fitType
  case {'greedy','greedy-local'} 
    res = loadFit(celldat,opts.fitType,[]);
    fit = res.fit_SUB;
    % get R2
    fit.r2_SUB = res.out_SUB.r2;
    fit.r2_LN = res.out_LN.r2;
    % get R2 on maximally differentiating stimuli
    celldat.loadType = 2;
    celldat.tempFilter = 1;
    celldat.zscore = 1;
    celldat.trainFrac = 0.8;
    celldat.subDivide = 10;
    [~,test,~,~] = loadCellData(dat,celldat,0);
    diffVals = (res.out_SUB.Z_t - res.out_LN.Z_t).^2;
    inds = diffVals > prctile(diffVals,80) & diffVals < prctile(diffVals,100);
    out_SUB_SEL = evalFit(test,res.fit_SUB,inds);
    out_LN_SEL =  evalFit(test,res.fit_LN,inds);
    fit.r2_SUB_SEL = out_SUB_SEL.r2;
    fit.r2_LN_SEL = out_LN_SEL.r2;
    
  case 'LN'
    res = loadFit(celldat,'greedy',[]);
    fit = res.fit_LN;
    fit.S = 1;
    fit.I_sc = ones(1,size(celldat.locs_c,1));
    fit.B_s = 0;
    fit.r2_LN = res.out_LN.r2;
end
