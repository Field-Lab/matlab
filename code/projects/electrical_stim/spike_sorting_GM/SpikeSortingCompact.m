function Output = SpikeSortingCompact(pathToAnalysisData,patternNo,neuronIds,pathToEi,varargin)
%SpikeSortingCompact does spike as shown in the examples (ExampleSeveralPatterns.m, ExampleSeveralNeurons.m)
%but allows to modify key parameters in the input. This function can be easily extended to allow more
%optional parameters (see below to see how optional parameters are modeled)
%No need to take care of misalignments in templates, those problems are solved
% authomatically using the AlingTemplates function
% input:    -pathToAnalysisData: path to ElecResp files and movie folders
%           -neuronIds: vector with Ids of neurons for which spike sorting is needed
%           -pathToEi: full path to the dataxxx.ei file
% optional: -recElecs: Electrodes that will be used for spike sorting. If no value is specified
%                      then electrodes from the electrode with the largest
%                      EI signal will be used for each neuron
%            -findAxon: logical, if zero no Axon breakpoints will be included (default =0)
%            -cleanData: if 1, the first trial of each movie will not be included. (defaul =1)
%            -collapseTrials: if 1, movies corresponding to same stimulus will be collapsed into one (default =1)
%            -Trange: Time range to do spike sorting, specified as a two dimensional vector (default = [1 40])
%            -degPolRule: polynomial rule for the Artifact covariates in initialization. If not specified, it will be the one coming from
%                         fillDefaultValues, that is, [1 4 15 50]
% Output:   - Output structure, (see  OutputDetails and documentation.rtf for details about this structure)
%If more optional parameters want to be used, notice it is important to figure out at which part of the code they have to be defined
% For example, degPolRule is defined after the default specification given in fillDefaultValues
% 
% 2015-07-13 LG edit: templates are generated directly from .ei files
% rather than taken from elecResp files (which do not necessarily exist).

cvx_solver Mosek
options = struct('recElecs',[],'findAxon',0,'cleanData',1,'collapseTrials',1,'Trange',[1 40],'degPolRule',[]);


nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('EXAMPLE needs propertyName/propertyValue pairs')
end


if(~isempty(varargin))
    optionNames = fieldnames(options);
    for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
        inpName = pair{1}; %# make case insensitive
        
        if any(strcmp(inpName,optionNames))
            options.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
end

input.tracesInfo.Trange        =    options.Trange;                                                                                       

if(isempty(options.recElecs))
    [templates, recElecs] = makeTemplatesFromEi(pathToEi,neuronIds); %This function will output one electrode per neuron
else
    recElecs = options.recElecs;
    [templates, recElecs] = makeTemplatesFromEi(pathToEi,neuronIds,recElecs);
end
fprintf('recording electrode for n%0.0f is %0.0f\n',neuronIds,recElecs);
templates = AlignTemplates(templates);
prefElectrodes = PreferredElectrodes(templates); % Looks for the template with the maximum signal for each neuron
input.neuronInfo.prefElectrodes =   prefElectrodes; %Define preferred electrodes
input.tracesInfo.recElecs      =    recElecs; %define recording electrodes
input.neuronInfo.neuronIds     =    neuronIds; %define neuron Ids
input.neuronInfo.templates     =    templates; %define templates


%Now Set defaults for how to load data. Should we collapse trials? should we
%clean the data first? Should we look for axonal activation (doesn't work for now)

%Parameters to find a Axonal activation authomatically. 
%if includeAxonBreakpoint=0 no axonal breakpoint will be included
input.params.load.findAxon.includeAxonBreakpoint  =   options.findAxon; % For now, set to zero because the axonal breakpoint solution is not working
input.params.load.findAxon.numberOmit             =   4;
input.params.load.findAxon.lMovingAverage         =   1;
input.params.load.findAxon.TrangeAxon             =  [11 40];
%If one, look for decreases in the minimum of the eigenvalues of the
%variance (space) of the energy plot. If equal 2 looks for increases
input.params.load.findAxon.typeEigenvalue         =   2;
% cleanData = 1 erases the first trial of each movie (for cases where some
% of the trials did not result in electrical stimulation as intended)
input.params.load.cleanData                       =   options.cleanData;
% if collapseTrialsSameCondition = 1 Collapse movies 2*j and 2*j-1
input.params.load.collapseTrialsSameCondition     =   options.collapseTrials;
 


%loads data from movie files
input = loadData(input,pathToAnalysisData,patternNo);

%fills with default values
input = fillDefaultValues(input);
if(~isempty(options.degPolRule))
    input.params.initial.degPolRule = options.degPolRule;
end
%any change in the default values after should be done at this point;
%input.params.initial.degPolRule{e}=[1 2 3 50]; example of a redefinition of a default value

%finds initial values for the artifact and rest of variables (via convex
%relaxation with cvx-Mosek)
initial = Initialize(input);
toc
disp('finished initialize');
%The Core of the algorithm: spike sorting via gibbs sampler + heuristics.
%solutions are stored in Gibbs Structure
[Gibbs, ~, ~, input, Log] = doSpikeSortingElectricalArtifact(input,initial);
toc
disp('finished doSpikeSortingElectricalArtifact');
%Creates the output structures
Output = OutputResults(input,Gibbs,Log);
   
    

