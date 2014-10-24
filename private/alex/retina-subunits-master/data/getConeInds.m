function celldat = getConeInds(dat,celldat)

sta = dat.coneWeights(:,celldat.cellInd);
if dat.stas.polarities{celldat.cellInd} == -1
	sta = -sta;
end
staMag = abs(sta);

weights = staMag/sum(staMag);
weightsTmp = weights;
weightsTmp(weights<0.5*max(weights)) = 0;
weightsTmp = weightsTmp/sum(weightsTmp);
center = weightsTmp'*dat.locations;
[distances diff] = getDistances(dat.locations,center);

keepInds = [];
switch celldat.threshType
  case 'distance'
    keepInds = distances < celldat.thresh;
  case 'peak'
    if dat.stas.polarities{celldat.cellInd} == -1
      keepInds = staMag > celldat.percent*max(staMag) & sta < 0 & distances < 40;
    end
    if dat.stas.polarities{celldat.cellInd} == 1
      keepInds = staMag > celldat.percent*max(staMag) & sta > 0 & distances < 40;
    end 
end

celldat.locs_c = dat.locations(keepInds,:);
celldat.coneIds = find(keepInds);
celldat.nConesFit = size(celldat.locs_c,1);
celldat.center = center;
celldat.sta = sta;
celldat.keepInds = keepInds;