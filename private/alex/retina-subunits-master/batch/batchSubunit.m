%function batchSubunit(computer,dataSet,cellType)
clear all
%%
computer = 'laptop';
dataSet = '2011-10-25-5/data001-0';
cellType = 'off midget'
dat = loadData(computer,dataSet);
plotting = 0;
%%
numCells = getCellTypeNum(dat.cellTypes,cellType);
for i=8:8;%numCells
    try
        fprintf('(batchSubunit) analyzing %s %g/%g\r',cellType,i,numCells);
        
        dat.cellType = 'off midget';
        dat.cellNum = i;
        dat.trainFrac = 0.8;
        dat.subDivide = 10;
        dat.thresh = 15;
        dat.percent = 0.2;
        dat.threshType = 'peak';
        dat.zscore = 1;
        dat.tempFilter = 1;
        [train test dat] = loadCellData(dat,1,0);
        
        doPairwise(dat);
        doSTC(dat,train,test,1);
        doSubunit(dat,train,test,'singletons',[],0);
        init_orig = fitGreedy(dat,train,test,1);
        doSubunit(dat,train,test,'greedy',init_orig,0);
        
        if plotting
            res = loadFit(dat,'greedy',20);
            res_sing = loadFit(dat,'singletons',50);
            fits.fit_SUB = res.fit_SUB;
            fits.fit_LN = res.fit_LN;
            fits.fit_SING = res.fit_SUB;
            plotModel2d(fits,dat,test);
        end
        
        close all
    catch
        fprintf('(batchSubunit) skipping %s %g/%g\r',cellType,i,numCells);
    end
end