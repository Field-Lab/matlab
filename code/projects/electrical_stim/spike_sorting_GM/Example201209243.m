cvx_solver Mosek
%pathToAnalysisData = './dataExample/';
% pathToAnalysisData = '/Users/gomena/Research/EJBigData/EJ-2014-11-05-Processed/data005/';
pathToAnalysisData = '/Volumes/Analysis/2012-09-24-3/data003/';
tic
neuronIds = [2318 2448 2209];
recElecs = [156 149 155];
patternNo = 156;
%recElecs  = 156;


%Default minimum and maximum recording times. If the first part of the
%traces are not relevant, it could be worth trying make Tmin>1
Tmin     = 1;
Tmax     = 40;

% Load Templates from ElecResp. If recElecs are not supplied, then will be
% found from the 'goodElecs' field.
[templates recElecs] = makeTemplatesFromElecResp(pathToAnalysisData,patternNo,neuronIds,recElecs,11);
%translate templates if by some reason there are weird template offsets.
%templates = translateTemplate(templates,-10,1,1);

%for each neuron, find the 'preFerred' electrode
prefElectrodes = PreferredElectrodes(templates);
input.neuronInfo.prefElectrodes =   prefElectrodes;

%Parameters to find a Axonal activation authomatically. 
%if includeAxonBreakpoint=0 no axonal breakpoint will be included
input.params.load.findAxon.includeAxonBreakpoint  =   0;
input.params.load.findAxon.numberOmit             =   4;
input.params.load.findAxon.lMovingAverage         =   1;
input.params.load.findAxon.TrangeAxon             =  [11 40];
%If one, look for decreases in the minimum of the eigenvalues of the
%variance (space) of the energy plot. If equal 2 looks for increases
input.params.load.findAxon.typeEigenvalue         =   2;
%cleanData = 1 erases the first trial of each movie
input.params.load.cleanData                       =   1;
%if collapseTrialsSameCondition = 1 Collapse movies 2*j and 2*j-1
input.params.load.collapseTrialsSameCondition     =   0;

input.tracesInfo.Trange        =    [Tmin Tmax];
input.tracesInfo.recElecs      =    recElecs;
input.neuronInfo.neuronIds     =    neuronIds;
input.neuronInfo.templates     =    templates;

input = loadData(input,pathToAnalysisData,patternNo);

input = filldefaultvalues201209243data3(input);

initial = Initialize(input);
toc
tic

[Gibbs GibbsNoDelete initial input Log] = doSpikeSortingElectricalArtifact(input,initial);
toc