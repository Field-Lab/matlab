function doPairwisePlot(celldat,dat)

% given a dat structure, grab the linearity matrix
% and save out a plot

if strcmp(celldat.computer,'salk');
    withNum = 1;
else
    withNum = 0;
end

v = loadVmat(celldat,dat);
plotPairwiseCell(celldat,v,withNum);
set(gcf, 'PaperPositionMode', 'auto');
if ~isdir(strcat(celldat.figurePath,'/pairwise/'))
    mkdir(strcat(celldat.figurePath,'/pairwise/'));
end
print('-dpdf',fullfile(celldat.figurePath,'pairwise',sprintf('%g-pairwise.pdf',celldat.rgcId)));