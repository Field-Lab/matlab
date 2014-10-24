function batchAnal(computer,dataSet,cellType,analType,doAnal,doPlot,percent,cellIds)

  dat = loadData(computer,dataSet);

  if exist('cellIds','var')
    loadType = 2;
    numCells = length(cellIds);
  else
    loadType = 1;
    numCells = getCellTypeNum(dat.cellTypes,cellType);
    cellIds = 1:numCells;
  end

  for i=1:length(cellIds)
    fprintf('(batchAnal) analyzing %s %g/%g from %s\n',cellType,find(cellIds(i)==cellIds),numCells,dat.dataSet);
    
    celldat = [];
    celldat.cellType = cellType;
    celldat.cellNum = cellIds(i);
    celldat.rgcId = cellIds(i);
    celldat.loadType = loadType;  
    celldat.percent = 0.2;
    celldat = getDefaultOpts(celldat);
    
    [train test celldat] = loadCellData(dat,celldat,0);
    if ~isempty(train) && ~isempty(test)
      if size(train.X_ct,1) > 1 && size(train.X_ct,1) < 100
        
        if sum(strcmp(analType,'pairwise'))
          if doAnal
            doPairwise(celldat,train,test);
          end
          if doPlot
            doPairwisePlot(celldat,dat);
          end
        end
        
        if sum(strcmp(analType,'stc'))
          if doAnal
            doSTC(celldat,train,test);
          end
          if doPlot
            doSTCplot(celldat,train,test);
          end
        end
        
        if sum(strcmp(analType,{'subunit','subunit-local','subunit-exhaustive'}))

          analMode = [];
          switch analType
          case 'subunit'
            analMode = 'greedy';
          case 'subunit-local'
            analMode = 'greedy-local';
          case 'subunit-exhaustive';
            analModel = 'exhaustive';
          end
          
          if percent ~= 0.2
            celldat.percent = percent;
            [train test celldat] = loadCellData(dat,celldat,0);
          end
          
          if ~isempty(train) && ~isempty(test)
            if doAnal
              res = loadFit(celldat,'LN');
              if ~isempty(res)
              if res.out_LN.r2 > 0.05 % r2 check for skipping subunit fit
                if size(train.X_ct,1) > 1 && size(train.X_ct,1) < 100
                  saveName = sprintf('%g-subunit-%s',celldat.rgcId,analMode);
                  saveName = strcat(saveName,'-',num2str(celldat.percent));
                  if ~isfile(fullfile(celldat.analPath,'subunit',strcat(saveName,'-fit.mat')));
                    doSubunit(celldat,train,test,'singletons',[]);
                    init_orig = [];
                    switch analMode
                    case 'greedy'
                      init_orig = fitGreedy(celldat,train,test,0);
                    case 'greedy-local'
                      init_orig = fitGreedyLocal(celldat,train,test,0);
                    case 'exhausitve'
                      init_orig = fitExhaustive(celldat,train,test,0);
                    end
                    doSubunit(celldat,train,test,analMode,init_orig);
                  else
                    fprintf('(batchAnal) fits exist for %g/%g\n',i,numCells);
                  end
                else
                  fprintf('(batchAnal) not enough cones, skipping %s %g/%g\n',cellType,i,numCells);
                end
              else
                fprintf('(batchAnal) r2 too low, skipping %s %g/%g\n',cellType,i,numCells);
              end
            else
              fprintf('(batchAnal) no LN fit, skipping %s %g/%g\n',cellType,i,numCells);
            end
          else 
            fprintf('(batchAnal) not enough cones, skipping %s %g/%g\n',cellType,i,numCells);
          end
        end
        
        if doPlot
          res = loadFit(celldat,analMode,[]);
          if ~isempty(res)
            saveName = sprintf('%g-subunit-%s',celldat.rgcId,analMode);
            saveName = strcat(saveName,'-',num2str(celldat.percent));
            savePath = sprintf('subunit-%s-%s-%s',analMode,num2str(celldat.percent),celldat.cellType);
            doSubunitPlot(celldat,dat,train,test,analMode);
          end
        end
        
      end
    else
      fprintf('(batchAnal) not enough cones, skipping %s %g/%g\n',cellType,i,numCells);
    end
  end
close all;
end