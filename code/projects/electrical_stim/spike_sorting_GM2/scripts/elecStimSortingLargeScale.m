function elecStimSortingLargeScale(pathToPreparation,dataFolders,pathSave,pathToEi,patternNos,neuronIds,templates,params)
% inputs:     pathToAnalysisData: '/Volumes/Analysis/...'choose the electrical stim data folder (001,002, etc) that contains data organized by pattern
%             pathToEi: a string, '/full/path/to/ei/file.ei'
%             patternNos
%             neuronIds


% codebase_path = matlab_code_path; 

% Set optional arguments. 


     for p = 1:length(patternNos)
        patternNo=patternNos(p);
        
          [elecRespAuto]=DoSpikeSortingLargeScaleNOEi(pathToPreparation,pathToEi,patternNo,neuronIds,params,templates,dataFolders);
       
   fname = fullfile(pathSave,['elecRespAuto_n' ...
        num2str(elecRespAuto.neuronInfo.neuronIds') '_p' ...
        num2str(elecRespAuto.stimInfo.patternNo) '.mat']); 
    save(fname,'elecRespAuto'); 
    disp(['done analyzing ' fname]); 
end




