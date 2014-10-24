function batchSimulate(computer,dataSet,cellType,percent,niters)

opts.durs = [5400, 10800, 21600]; % 7.5 min, 15 min, 30 min
%durs = [2880, 5760, 11520, 23040]; % 4 min, 8 min, 16 min, 32 min
opts.cellIds = [677 302 481 1066 1247 7252];

dat = loadData(computer,dataSet);

for icell=1:length(opts.cellIds)
	celldat = [];
	cellType = 'on midget';
	celldat.rgcId = opts.cellIds(icell);
	celldat.loadType = 2;  
	celldat.percent = percent;
	celldat = getDefaultOpts(celldat);
	[train test celldat] = loadCellData(dat,celldat,0);
	res = loadFit(celldat,'greedy',[]);
	origfit = res.fit_SUB;

	for idur=1:length(opts.durs)
		fprintf('(batchSimulate) simulating cell %g/%g, duration %g\n',icell,length(opts.cellIds),opts.durs(idur));
		for iiter=1:niters
			fprintf('(batchSimulate) simulation %g/%g\n',iiter,niters);
			[simdat simtrain simtest] = getSimulation(celldat,train,test,origfit,opts.durs(idur),1000);
			simdat.simopts.iter = iiter;

			% fit to the simulated spikes (greedy)
			simdat.simopts.mode = '';
			analMode = 'greedy';
			init_orig = fitGreedy(simdat,simtrain,simtest,0);
			doSubunit(simdat,simtrain,simtest,analMode,init_orig);

			% fit to the simulated spikes (correct structure)
			simdat.simopts.mode = 'fixed';
			analMode = 'greedy';
			init_orig = origfit.I_sc;
			doSubunit(simdat,simtrain,simtest,analMode,init_orig);
		end
	end
end
