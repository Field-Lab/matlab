function batchVmat(computer,dataset,strt)

dat = loadData(computer,dataset);

if ~exist('strt','var') || isempty(strt)
    strt = 1;
end

nCells = length(dat.rgcIds);
for iCell=strt:nCells
    if ~isempty(dat.stas.time_courses{iCell})
        dat.cellNum = dat.rgcIds(iCell);
        dat.trainFrac = 0.8;
        dat.subDivide = 10;
        dat.thresh = 2;
        dat.percent = 0.2;
        dat.threshType = 'peak';
        dat.zscore = 1;
        dat.tempFilter = 1;
        dat.crossVal = 0;
        [train test dat] = loadCellData(dat,2,0);
        %%
        if ~isempty(train) && ~isempty(test)
            if size(train.X_ct,1) > 0 && size(train.X_ct,1) < 100
                if strcmpi(dat.cellType,'on midget') || ...
                        strcmpi(dat.cellType,'off midget') || ...
                        strcmpi(dat.cellType,'on parasol') || ...
                        strcmpi(dat.cellType,'off parasol') || ...
                        strcmpi(dat.cellType,'sbc')
                    
                    doPairwise(dat,train,test);
                end
            end
        end
    end
end