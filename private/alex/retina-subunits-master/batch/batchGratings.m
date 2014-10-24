function batchGratings(computer,dataSet,cellType,percent,datarun)

dat = loadData(computer,dataSet);
numCells = getCellTypeNum(dat.cellTypes,cellType);
cellIds = 1:numCells;

expvar_LN = [];
expvar_SUB = [];
expvar_LN_per = [];
expvar_SUB_per = [];

for icell=1:numCells
	icell
	celldat = [];
	celldat.cellNum = icell;
	celldat.loadType = 1;  
	celldat.cellType = cellType;
	celldat.percent = percent;
	celldat = getDefaultOpts(celldat);
	[train test celldat] = loadCellData(dat,celldat,0);
	res = loadFit(celldat,'greedy',[]);
	g = getGratingResp2(datarun{2},celldat);

	if ~isempty(g) && ~isempty(res)
		
		fit_LN = res.fit_LN;
		fit_SUB = res.fit_SUB;

		opts.pers = [4:8]; % was using 4:8
		gfit = getGratingFit(g,fit_LN,opts);
		expvar_LN = [expvar_LN; gfit.cr];

		gfit = getGratingFit(g,fit_SUB,opts);
		expvar_SUB = [expvar_SUB; gfit.cr];

		expvar_LN_per_tmp = [];
		expvar_SUB_per_tmp = [];
		for iper=4:8
			opts.pers = iper;
			gfit = getGratingFit(g,fit_LN,opts);
			expvar_LN_per_tmp = [expvar_LN_per_tmp, gfit.cr];
			gfit = getGratingFit(g,fit_SUB,opts);
			expvar_SUB_per_tmp = [expvar_SUB_per_tmp, gfit.cr];
		end
		expvar_LN_per = [expvar_LN_per; expvar_LN_per_tmp];
		expvar_SUB_per = [expvar_SUB_per; expvar_SUB_per_tmp];
	end
end

length(expvar_LN)

clf
set(gcf,'Position',[0 0 950 400]);

subplot(1,3,1);
bins = linspace(0,1,20);
h = hist(expvar_SUB,bins);
b = bar(bins,h/sum(h));
xlim([-0.05,1.05]);
axis square;
box off; set(gca,'TickDir','out');
shading flat
set(gca,'FontSize',18);
ylabel('Fraction of cells');
xlabel('Explained variance');
set(b,'FaceColor',[0.5 0.5 0.5]);
set(b,'EdgeColor','none');
title(sprintf('n = %g',length(expvar_SUB)));

subplot(1,3,2);
plot(expvar_LN,expvar_SUB,'k.','MarkerSize',20);
axis square;
xlim([0 1]);
ylim([0 1]);
set(gca,'XTick',[0 0.5 1]);
set(gca,'YTick',[0 0.5 1]);
set(gca,'FontSize',18);
xlabel('Adjusted R2 (LN)');
ylabel('Adjusted R2 (Subunit)');
box off; set(gca,'TickDir','out');
hold on;
plot([0 1],[0 1],'--');

%badinds = zeros(1,36);
%badinds([13,23,24,25]) = 1;

subplot(1,3,3);
LN_mn = nanmean(expvar_LN_per);
LN_sd = bsxfun(@rdivide,nanstd(expvar_LN_per),sqrt(sum(~isnan(expvar_LN_per))));
SUB_mn = nanmean(expvar_SUB_per);
SUB_sd = bsxfun(@rdivide,nanstd(expvar_SUB_per),sqrt(sum(~isnan(expvar_SUB_per))));
plotErrorBars(g.perVals(4:8),LN_mn(1:end),LN_sd(1:end),[0.5 0.5 0.5],'b',0,1);
hold on
h = plot(g.perVals(4:8),LN_mn(1:end),'o');
set(h,'MarkerEdgeColor',[1 1 1]);
set(h,'MarkerFaceColor',[0 0 1]);
set(h,'MarkerSize',12);
plotErrorBars(g.perVals(4:8),SUB_mn(1:end),SUB_sd(1:end),[0.5 0.5 0.5],'g',0,1);
hold on
h = plot(g.perVals(4:8),SUB_mn(1:end),'o');
set(h,'MarkerEdgeColor',[1 1 1]);
set(h,'MarkerFaceColor',[0 1 0]);
set(h,'MarkerSize',12);
set(gca,'FontSize',18);
set(gca,'XScale','log');
xlim([3 25]);
ylim([0 1]);
ylabel('Adjusted R2');
xlabel('Spatial frequency (cycles/deg)');
box off
set(gca,'TickDir','out');
axis square;
xlim([4 50]);

keyboard

% subplot(2,1,2);
% improvement = (expvar_SUB_per)./(expvar_LN_per);
% improvement_mn = nanmean(improvement);
% improvement_sd = nanstd(improvement)/sqrt(length(expvar_SUB));
% plotErrorBars(g.perVals(7:g.nPer-1),improvement_mn(1:end-1),improvement_sd(1:end-1),[0.5 0.5 0.5],'k',0,1);
% hold on
% h = plot(g.perVals(7:g.nPer-1),improvement_mn(1:end-1),'o');
% set(h,'Color','k');
% set(h,'MarkerEdgeColor',[1 1 1]);
% set(h,'MarkerFaceColor',[0 0 0]);
% set(h,'MarkerSize',12);
% set(gca,'FontSize',14);
% set(gca,'XScale','log');
% xlim([3 25]);
% ylim([0.7 2.5]);