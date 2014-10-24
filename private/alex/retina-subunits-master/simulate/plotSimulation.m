function plotSimulation(dat,rgcId,percent,niters)

	clf

opts.durs = [5400, 10800, 21600]; % 7.5 min, 15 min, 30 min
celldat = [];
celldat.rgcId = rgcId;
celldat.loadType = 2;  
celldat.percent = percent;
celldat = getDefaultOpts(celldat);
[train test celldat] = loadCellData(dat,celldat,0);
res = loadFit(celldat,'greedy',[]);
origfit = res.fit_SUB;

coneSz = 1.4;
subSz = 2;
coneEdge = 1.25;

sanesubplot(length(opts.durs),niters+3,[1 1],0.02);
axis equal
xlim(celldat.xRange);
ylim(celldat.yRange);
set(gca,'YDir','reverse');
axis off
%rect1 = rectangle('Position',[celldat.xRange(1) celldat.yRange(1) sum(celldat.xRange) sum(celldat.yRange)]);
%set(rect1,'FaceColor',[0.5 0.5 0.5]);
%set(rect1,'EdgeColor','none');
plotSubunitsLinesOvals(origfit,celldat,coneSz,subSz,coneEdge);

for idur=1:length(opts.durs)

	for iiter=1:niters
		[simdat, ~, ~] = getSimulation(celldat,train,test,origfit,opts.durs(idur),1000);
		simdat.simopts.iter = iiter;
		simdat.simopts.mode = '';
		res = loadFit(simdat,'greedy',[]);
		thisfit = res.fit_SUB;

		% plot the subunits for each iter as columns
		a = sanesubplot(length(opts.durs),niters+3,[idur iiter+1],0.02);
		axis equal
		xlim(celldat.xRange);
		ylim(celldat.yRange);
		set(gca,'YDir','reverse');
		axis off
		%rect1 = rectangle('Position',[celldat.xRange(1) celldat.yRange(1) sum(celldat.xRange) sum(celldat.yRange)]);
		%set(rect1,'FaceColor',[0.5 0.5 0.5]);
		%set(rect1,'EdgeColor','none');
		plotSubunitsLinesOvals(thisfit,celldat,coneSz,subSz,coneEdge);

		simdat.simopts.mode = 'fixed';
		res = loadFit(simdat,'greedy',[]);
		thisfit = res.fit_SUB;

		if iiter == 1
			h1 = sanesubplot(length(opts.durs),niters+3,[idur niters+2],0.02);
		end

		% should only be for subunits with more than one cone...
		plot(h1,origfit.A_sc(logical(origfit.I_sc)),thisfit.A_sc(logical(origfit.I_sc)),'b.');
		hold on
		plot(h1,origfit.B_s,thisfit.B_s,'r.');
		set(h1,'Box','off');
		axis(h1,'equal');
		xlim(h1,[-0.25 1.25]);
		ylim(h1,[-0.25 1.25]);
		plot(h1,[-1 2],[-1 2],'k--');
		set(h1,'XTick',[0 0.5 1]);
		set(h1,'YTick',[0 0.5 1]);
		set(h1,'TickDir','out');

		if iiter == 1
			h2 = sanesubplot(length(opts.durs),niters+3,[idur niters+3],0.01);
		end

		axes(h2)
		if iiter == 1
		 	plotNonLin(origfit.f,'r','line');
		end
		plotNonLin(thisfit.f,'k','line');

	end
end
