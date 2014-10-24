function out = loadVmat(celldat,dat)

savePath = strcat(celldat.analPath,'/','vMats','/');
fileName = sprintf('vmat-%g',celldat.rgcId);

if isfile(fullfile(savePath, [fileName,'.mat']));
	
	load(fullfile(savePath,fileName));
	out = v.mat;

	if celldat.percent ~= v.percent
	    orig_inds = celldat.coneIds;
	    celldat2 = celldat;
	    celldat2.percent = v.percent;
	    [~, ~, celldat2] = loadCellData(dat,celldat2,0);
	    [~,~,ib] = intersect(orig_inds,celldat2.coneIds);
	    out = out(ib,ib);
	end
else
	out = [];
end
