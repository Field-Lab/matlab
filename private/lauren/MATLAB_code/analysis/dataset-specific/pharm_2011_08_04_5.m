
clear all

cellIDsVis =   [291      571];
cellIDsStim =  [291      571];
stimElecs =    [17       43];
cellTypes =    {'offP', 'onP'};

phaseWidth = 100;

pathToElecRespFilesBefore = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data009/';
pathToElecRespFilesAfter  = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data011/';

pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data010/data010.neurons';

phaseWidth = 100;


%% plot summary figure for each cell

for ii = 1:length(cellIDsStim)
    pathToElecResp.before = [pathToElecRespFilesBefore 'elecResp_n'...
        num2str(cellIDsStim(ii)) '_p' num2str(stimElecs(ii)) '_w' num2str(phaseWidth)];
    pathToElecResp.after  = [pathToElecRespFilesAfter  'elecResp_n'...
        num2str(cellIDsStim(ii)) '_p' num2str(stimElecs(ii)) '_w' num2str(phaseWidth)];
    
    pharm_results_plotter(cellIDsVis(ii), cellTypes{ii}, pathToElecResp,...
        pathToNeuronsFile)
end

