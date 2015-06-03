%path('/Users/gomena/Research/GIT/spike-sorting-electrical-artifact')
addpath('./Core/')
addpath('./Utils/')
cvx_solver Mosek
%pathToAnalysisData = './dataExample/';
pathToAnalysisData = '/Users/gomena/Research/EJ-2014-11-05-Processed/data005/';

sampRate           = 20000;

neuronIds = [2418 2434 2913 2763];
patternNo = 1443;
recElecs  = [170 233 194 185];

Tmin     = 1;
Tmax     = 40;
Tfindmin = 5;
Tfindmax = 35;
T        = Tmax-Tmin+1;
Tdivide  = [0 floor(T/2) T];

templates = makeTemplatesFromElecResp(pathToAnalysisData,patternNo,neuronIds,recElecs);
templates = translateTemplate(templates,6,2,2);

movieNos  = findMovieNos(pathToAnalysisData,patternNo);

%Translate templates if necessary

   

for m = 1:size(movieNos,2)
    
    dataTraces = NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
            movieNos(m), 99999);
      
    [amps channelsWithStim stimAmpVectors channelsConnected elecCurrentStep currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
     
     listAmps(m,:)=amps';
     listStimElecs(m,:)=channelsWithStim;
     listCurrentRangesUsed(m,:)=currentRangesUsed;
     
    
     for e = 1:length(recElecs)
            
            data{m,e}=squeeze(dataTraces(:,recElecs(e),Tmin:Tmax));
     
     end
end


[data listAmps listStimElecs  listCurrentRangesUsed] = collapseTrialsSameCondition(data,listAmps,listStimElecs,listCurrentRangesUsed);

[data] = cleanTrials(data);
dataVecJ  =  collapseDataE(data);
%dont forget there are patholicical cases, with collapsing + breakpoints
breakStimElecs = findBreakStimElecs(listCurrentRangesUsed);
%explore what happens with breakpoints in other datsets
[breakRecElecs]  = findBreakRecElecs(breakStimElecs,recElecs,listStimElecs);


input.tracesInfo.data          =    data;
input.tracesInfo.dataVecJ      =    dataVecJ;
input.tracesInfo.Trange        =    [Tmin Tmax];
input.tracesInfo.recElecs      =    recElecs;
input.tracesInfo.breakRecElecs =    breakRecElecs;
input.stimInfo.listAmps        =    listAmps;
input.stimInfo.breakStimElecs  =    breakStimElecs;
input.stimInfo.listStimElecs   =    listStimElecs;
input.stimInfo.patternNo       =    patternNo;
input.stimInfo.movieNos        =    movieNos;
input.neuronInfo.neuronIds     =    neuronIds;
input.neuronInfo.templates     =    templates;
input.params.sampRate          =    sampRate;
input.params.path              =    pathToAnalysisData;
input.params.TfindRange        =    [Tfindmin Tfindmax];
input.params.TfindRangeRel     =    input.params.TfindRange - input.tracesInfo.Trange(1) +1;
input.params.Tdivide           =    Tdivide;

input = fillDefaultValues(input);

initial = Initialize(input);

[Gibbs initial input] = doSpikeSortingElectricalArtifact(input,initial);
