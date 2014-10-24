dat = loadData('laptop-2','2010-03-05-2/data013');
celldat = [];
cellType = 'on midget';
celldat.rgcId = 677;
celldat.loadType = 2;  
celldat.percent = 0.33;
celldat = getDefaultOpts(celldat);
[train test celldat] = loadCellData(dat,celldat,0);
res = loadFit(celldat,'greedy',[]);
fit = res.fit_SUB;

[simdat simtrain simtest] = getSimulation(celldat,train,test,origfit,1000,1000);

% fit to the simulated spikes (greedy)
simdat.simopts.mode = '';
analMode = 'greedy';
init_orig = fitGreedy(simdat,simtrain,simtest,1);
doSubunit(simdat,simtrain,simtest,analMode,init_orig);

% fit to the simulated spikes (correct structure)
simdat.simopts.mode = 'fixed';
analMode = 'greedy';
init_orig = origfit.I_sc;
doSubunit(simdat,simtrain,simtest,analMode,init_orig);

% do pairwise on the simulation
%doPairwise(simdat,sim.train,sim.test);
% fit to the simulated spikes (singleton)
%doSubunit(simdat,sim.train,sim.test,'singletons',[]);