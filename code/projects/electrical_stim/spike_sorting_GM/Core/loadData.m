function input = loadData(input,pathToAnalysisData,patternNo)

input.names.path               =    pathToAnalysisData;
input.stimInfo.patternNo       =    patternNo;

recElecs                    = input.tracesInfo.recElecs;
Trange                      = input.tracesInfo.Trange;
neuronIds                   = input.neuronInfo.neuronIds;
TrangeAxon                  = input.params.load.findAxon.TrangeAxon;
includeAxonBreakpoint       = input.params.load.findAxon.includeAxonBreakpoint;
cleanData                   = input.params.load.cleanData;
collapseTrialsSameCondition = input.params.load.collapseTrialsSameCondition;

movieNos  = findMovieNos(pathToAnalysisData,patternNo);
movieNos  = sort(movieNos);



dataTraces    = NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
                movieNos(1), 99999);
firstArtifact = mean(dataTraces,1);

for m = 1:length(movieNos)
    
    dataTraces = NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
            movieNos(m), 99999);
      
    
    subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);

    PowerMatrix = 0;
    meanT = mean(mean(dataTraces(:,:,TrangeAxon(1):TrangeAxon(2))-subtractionMatrix(:,:,TrangeAxon(1):TrangeAxon(2)),3),1);
    
    [~,EIm_viewT]   = ei2matrix(meanT);
    for t = TrangeAxon(1):TrangeAxon(2)
      
    	meanData     = mean(dataTraces(:,:,t)-subtractionMatrix(:,:,t),1);
        [~,EIm_view] = ei2matrix(meanData);
        PowerMatrix  = PowerMatrix+(EIm_view-EIm_viewT).^2;
       
    end
        
    PowerMatrix  = PowerMatrix; 
    
    Powers{m}    = PowerMatrix;
    
    
    for e = 1:length(recElecs)
            
        data{m,e}=squeeze(dataTraces(:,recElecs(e),Trange(1):Trange(2)));
     
    end    
    [amps channelsWithStim stimAmpVectors channelsConnected elecCurrentStep currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
     
     listAmps(m,:)              = amps';
     listStimElecs(m,:)         = channelsWithStim;
     listCurrentRangesUsed(m,:) = currentRangesUsed;
     

end

input.tracesInfo.Powers = Powers;
breakStimElecs          = findBreakStimElecs(listCurrentRangesUsed);
[breakRecElecs]         = findBreakRecElecs(breakStimElecs,recElecs,listStimElecs);


if includeAxonBreakpoint
    [input AxonBr]      = findAxonalBreakpoint(input);
else
    AxonBr      = [];
end

for e=1:length(recElecs)
        
    breakAxon{e}    = AxonBr; 
    breakPoints{e}  = unique(sort([breakRecElecs{e} breakAxon{e}]));

end


if cleanData
[data] = cleanTrials(data);
end

if collapseTrialsSameCondition
    [data listAmps listStimElecs  listCurrentRangesUsed] = collapseTrialsSameCondition(data,listAmps,listStimElecs,listCurrentRangesUsed);
end

dataVecJ  =  collapseDataE(data);

input.tracesInfo.data          =    data;
input.tracesInfo.dataVecJ      =    dataVecJ;
input.tracesInfo.breakAxon     =    breakAxon;
input.tracesInfo.breakRecElecs =    breakRecElecs;
input.tracesInfo.breakPoints   =    breakPoints;
input.stimInfo.listAmps        =    listAmps;
input.stimInfo.breakStimElecs  =    breakStimElecs;
input.stimInfo.listStimElecs   =    listStimElecs;

