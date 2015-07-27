function input = loadData(input,pathToAnalysisData,patternNo)
% loadData loads data from a movies adds manually an Axonal Bundle Breakpoint to the input structure, defined as the
% last condition j for which there is no observed traveling wave following stimulation.
% inputs:   -pathToAnalysisData: path
%           -patternNo: pattern to be analyzed
%           -input: input structure with a few values values (input.tracesInfo.Trange, input.params.findAxon. (see examples.m
%                   input.params.load.cleanData, input.params.load.collapseTrialsSameCondition )
%           
% output:  -input The same input, with data loaded (new information in fields
%           input.tracesInfo, input.stimInfo). Traces will be available in
%           input.tracesInfo.data{e,j}
% usage:  usit it in the same order as in examples (some fields in input have to be filled)
%
% Gonzalo Mena 6/2015 


input.names.path               =    pathToAnalysisData;
input.stimInfo.patternNo       =    patternNo;

recElecs                    = input.tracesInfo.recElecs;
Trange                      = input.tracesInfo.Trange;
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
    
    if includeAxonBreakpoint
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
        
    end
    for e = 1:length(recElecs)
        
        data{m,e}=squeeze(dataTraces(:,recElecs(e),Trange(1):Trange(2)));
        
    end
    [amps, channelsWithStim, ~, ~, ~, currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
    
    listAmps(m,:)              = amps';
    listStimElecs(m,:)         = channelsWithStim;
    listCurrentRangesUsed(m,:) = currentRangesUsed;
    
    
end
if includeAxonBreakpoint
    input.tracesInfo.Powers = Powers;
end

if cleanData
    [data] = cleanTrials(data);
end

if collapseTrialsSameCondition

    index1=[1:size(data,1)];
    cont=1;
    
    while(~isempty(index1))
        
        amps = listAmps(index1(1),:);
        
        
        indequal = strmatch(amps,listAmps)';
        
        
        for e=1:size(data,2)
            
            dataaux{cont,e} = [];
            
            for l=indequal
                
                dataaux{cont,e} = [dataaux{cont,e}; data{l,e}];
                
            end
            
        end
        
        listAmpsAux(cont,:)                =  listAmps(index1(1),:);
        listStimElecsAux(cont,:)            =  listStimElecs(index1(1),:);
        listCurrentRangesUsedAux(cont,:)   =  listCurrentRangesUsed(index1(1),:);
        
        index1=setdiff(index1,indequal);
        cont=cont+1;
        
    end

data                  = dataaux;
listCurrentRangesUsed = listCurrentRangesUsedAux;
listStimElecs          = listStimElecsAux;
listAmps               = listAmpsAux;


end



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


dataVecJ  =  collapseDataE(data);

input.tracesInfo.data          =    data;
input.tracesInfo.dataVecJ      =    dataVecJ;
input.tracesInfo.breakAxon     =    breakAxon;
input.tracesInfo.breakRecElecs =    breakRecElecs;
input.tracesInfo.breakPoints   =    breakPoints;
input.stimInfo.listAmps        =    listAmps;
input.stimInfo.breakStimElecs  =    breakStimElecs;
input.stimInfo.listStimElecs   =    listStimElecs;

