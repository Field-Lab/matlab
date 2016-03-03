
%% 
% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('../../../../matlab_scratch/matlab/'));

%% Load cell IDs
WN_datafile = '2015-11-09-3/data000/data000';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun = load_sta(datarun);


plot_rf_summaries(datarun,'ON Parasol');

stas=[];cellID_list=[];
for cellType=1:2
for icell=1:length(datarun.cell_types{cellType}.cell_ids)
    cellID = datarun.cell_types{cellType}.cell_ids(icell);
    ssta = datarun.stas.stas{datarun.cell_ids==cellID};
    ssta = mean(ssta(:,:,:,25),3);
    stas=[stas;ssta(:)'];
    cellID_list = [cellID_list;cellID];
end
end
%% load elecResp and make sigmoid
nCells = length(cellID_list);
sigmoid_params1 = zeros(nCells,512);
sigmoid_params2 = zeros(nCells,512);
icell=0;
for cellID = cellID_list'
    icell = icell+1
    try
    files = ls(sprintf('/Volumes/Analysis/2015-11-09-3/data001-data002-autosort/elecRespAuto_n%d*',cellID));
    files = strsplit(files,'\n');
    
    for ifile=1:length(files)-1
        
        elecResp = load(files{ifile});
        electrode = elecResp.elecRespAuto.stimInfo.listStimElecs(1);
        
        curve = (elecResp.elecRespAuto.LogisticReg);
        
        weight = 1:length(curve);
        [logitCoef,dev] = glmfit(weight,[curve' ones(length(curve),1)],'binomial','logit');
        logitFit = glmval(logitCoef,weight,'logit');
        plot(weight,curve,'bs', weight,logitFit,'r-');
        
        sigmoid_params1(icell,electrode)=logitCoef(1);
        sigmoid_params2(icell,electrode)=logitCoef(2);
    end
    catch 
    end
end


save('/Volumes/Lab/Users/bhaishahster/2015-11-09-3_data000_2.mat','stas','cellID_list','sigmoid_params1','sigmoid_params2','-v7.3');
