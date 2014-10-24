computer = 'laptop';

dat = loadData(computer,'2012-09-21-2/data007'); % wn run
load('~/Dropbox/Projects/retinaSubunits/data/2012-09-21-2/data008/conepreprocess.mat'); % repeat
% good inds: 6 (some big effects), 7 (a few big effects), 103, 116 (small distributed effects), 121 (big effects)

%dat = loadData(computer,'2012-09-18-3/data003'); % wn run
%load('~/Dropbox/Projects/retinaSubunits/data/2012-09-18-3/data004/conepreprocess.mat'); % repeat

%dat = loadData(computer,'2012-08-21-0/data001'); % wn run
%load('~/Dropbox/Projects/retinaSubunits/data/2012-08-21-0/data002/conepreprocess.mat'); % repeat


%%
%'2012-08-21-0/data001' 767, 2072, 2971, 4561, 5567, 6949
%'2012-09-18-3/data003' 422, 799, 993, 1517, 1636, 2868, 3137, 3241 
%                       4351, 7488
%'2012-09-21-2/data007' 91, 181, 331, 332, 1186

%best = [2666, 3318, 7401];
%best = [5567, 5852, 6286, 6602, 6814, 6668, 6949];

cr_all = [];
cr_LN_all = [];
best_cr_all = [];

cellType = 'off midget';
numCells = getCellTypeNum(dat.cellTypes,cellType);
cellIds = 1:numCells;

for i = 1:numCells;

	celldat = [];
	celldat.cellType = cellType;
	celldat.cellNum = cellIds(i);
	celldat.rgcId = cellIds(i);
	celldat.loadType = 1;  
	celldat.percent = 0.33;
	celldat = getDefaultOpts(celldat);

	[train test celldat] = loadCellData(dat,celldat,0);
	re = getRepeat(celldat,datarun);
	
	if ~isempty(re)
		%plotRepeat(re,celldat);
		%i,
		%pause
	end

	if ~isempty(re)
		cr_all = [cr_all,re.cr_boot/re.reliability];
		cr_LN_all = [cr_LN_all,re.cr_LN_boot/re.reliability];
	end
	
end

length(cr_all)

clf
set(gcf,'Position',[0 0 950 400]);

subplot(1,3,1);
bins = linspace(0,1,18);
h = hist(cr_all,bins);
b = bar(bins,h/sum(h));
xlim([-0.05,1.05]);
axis square;
%drawVertLine(best_cr_all);
box off; set(gca,'TickDir','out');
shading flat
set(gca,'FontSize',18);
ylabel('Fraction of cells');
xlabel('Adjusted R2');
set(b,'FaceColor',[0.5 0.5 0.5]);
set(b,'EdgeColor','none');

subplot(1,3,2);
plot(cr_LN_all,cr_all,'k.','MarkerSize',20);
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
