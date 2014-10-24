% for the pixelated noise, each stimulus was 8 frames
% for the gratings, each stimulus was 1 frame
% so would need to bin at (8/120) to get into alignment (8 frames)
dat = loadData('laptop','2010-03-05-2/data013'); % wn run
load /Users/freemanj11/Dropbox/Projects/retinaSubunits/data/2010-03-05-2/data016/datarun.mat

% 5, 8, 26, 29* -> good for plotting
% 13, 23, 24, 25, 26

for i=[29]
cellType = 'off midget';
numCells = getCellTypeNum(dat.cellTypes,cellType);
i = i;
celldat = [];
celldat.cellType = cellType;
celldat.cellNum = i;
celldat.rgcId = 677;
celldat.loadType = 1;  
celldat.percent = 0.33;
celldat = getDefaultOpts(celldat);
[train test celldat] = loadCellData(dat,celldat,0);

disp(celldat.rgcId)

g = getGratingResp(datarun{2},celldat);
plotGratings(g,datarun{2},celldat,'cycle');


if ~isempty(g)
	res = loadFit(celldat,'greedy',[]);
	if ~isempty(res)
		fit = res.fit_SUB;
		opts.pers = [4:8]; % was using 4:8
		[gfit refit] = getGratingFit(g,fit,opts);
		plotGratings(g,datarun{2},celldat,'cycle',gfit);
	end
end

pause;
end