function out = batchRetinaSummary(dataSet,cellType,analType)

% out = batchRetinaSummary(dat)
%
% get summary statistics for a retina in a common format

fprintf('(batchRetinaSummary) loading %s\n',dataSet);

dat = loadData('laptop',dataSet);

killList = loadKillList(dataSet,cellType);

coneThresh = 0.33;
r2Thresh = 0.05;

numCells = getCellTypeNum(dat.cellTypes,cellType);

for i=1:numCells
  
  fprintf('(batchRetinaSummary) getting %s %g/%g\n',cellType,i,numCells);
  
  celldat = [];
  celldat.cellType = cellType;
  celldat.cellNum = i;
  celldat.rgcId = [];
  celldat.loadType = 1;
  celldat.trainFrac = 0.8;
  celldat.subDivide = 10;
  celldat.thresh = 15;
  celldat.percent = coneThresh;
  celldat.threshType = 'peak';
  celldat.zscore = 1;
  celldat.tempFilter = 1;
  [train, test, celldat] = loadCellData(dat,celldat,0);
  
  if ~isempty(celldat)
    if sum(celldat.rgcId==killList)
      fprintf('(batchRetinaSummary) killing cell %g \n',celldat.rgcId);
      celldat = [];
    end
  end


  % get subunit counts
  if sum(strcmp(analType,'counts'))
    if ~isempty(celldat)
      res = loadFit(celldat,'greedy',[]);
      if ~isempty(res)
        %if (res.out_SUB.r2 >= r2Thresh)
          subSz = sum(res.fit_SUB.I_sc,2);
          for j=1:10
            counts.subCount(i,j) = sum(subSz==j);
          end
          counts.density(i,1) = size(res.fit_SUB.I_sc,2); % number of cones
          counts.density(i,2) = size(res.fit_SUB.I_sc,1); % number of subunits
        %end
      end
    else
      counts.subCount(i,:) = NaN;
      counts.density(i,1) = NaN;
      counts.denisty(i,2) = NaN;
    end
  end
  
  if sum(strcmp(analType,'nonlin'))
    if ~isempty(celldat)
      res = loadFit(celldat,'greedy',[]);
      resSING = loadFit(celldat,'singletons',[]);
      if ~isempty(res)
        diffVals = (res.out_SUB.Z_t - res.out_LN.Z_t).^2;
        inds = diffVals > prctile(diffVals,80) & diffVals < prctile(diffVals,100);
        out_SUB_SEL = evalFit(test,res.fit_SUB,inds);
        out_LN_SEL =  evalFit(test,res.fit_LN,inds);
        nonlin.r2.ln(i) = res.out_LN.r2;
        nonlin.r2select.ln(i) = out_LN_SEL.r2;
        nonlin.r2.sub(i) = res.out_SUB.r2;
        nonlin.r2.sing(i) = resSING.out_SUB.r2;
        nonlin.r2select.sub(i) = out_SUB_SEL.r2;
        nonlin.f.x = linspace(-2.25,2.25,100);
        nonlin.g.x = linspace(-1,7,100);
        %if (res.out_SUB.r2 >= r2Thresh)
          nonlin.f.out(:,i) = evalNonLin(nonlin.f.x,res.fit_SUB.f);
          nonlin.g.out(:,i) = evalNonLin(nonlin.g.x,res.fit_SUB.g);
        %end
      end
    else
      nonlin.r2.ln(i) = NaN;
      nonlin.r2select.ln(i) = NaN;
      nonlin.r2.sub(i) = NaN;
      nonlin.r2.sing(i) = NaN;
      nonlin.r2select.sub(i) = NaN;
    end
  end

  % get subunit blur info
  if sum(strcmp(analType,'blury'))
    if ~isempty(celldat)
      res = loadFit(celldat,'greedy',[]);
      if ~isempty(res)
        if res.out_SUB.r2 >= r2Thresh;
          if size(celldat.rf,3) > 1
            celldat.rf = mean(celldat.rf,3);
          end
          
          celldat = fitGaussCones(celldat);
          dists = getDistances(celldat.locs_c,celldat.locs_c);
          sz_sd = sqrt(celldat.rfFit.sz);
          [X Y] = meshgrid(sz_sd);
          sz_pairs = X+Y;
          normdists = 2*dists./(sz_pairs);
          subness = res.fit_SUB.I_sc'*res.fit_SUB.I_sc;
          
          % normalize the assignment matrix
          foo = normMat(res.fit_SUB.I_sc,'L2',2);
          % only entries that equal 1 are singletons
          singInds = sum(foo)==1;
          nonSingInds = ~singInds;
          
          blury.normdists.subpairs{i} = normdists(triu(logical(subness),1));
          blury.normdists.nonsubpairs{i} = normdists(triu(logical(1-subness),1));
          blury.dists.subpairs{i} = dists(triu(logical(subness),1));
          blury.dists.nonsubpairs{i} = dists(triu(logical(1-subness),1));
          blury.sz.subpairs{i} = sz_sd(nonSingInds);
          blury.sz.nonsubpairs{i} = sz_sd(singInds);
        end
      end
    end
  end
  
end



if sum(strcmp(analType,'counts'))
  out.counts = counts;
end
if sum(strcmp(analType,'blury'))
  out.blury = blury;
end
if sum(strcmp(analType,'nonlin'))
  out.nonlin = nonlin;
end



