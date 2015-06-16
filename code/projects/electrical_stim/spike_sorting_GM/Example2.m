cvx_solver Mosek
%pathToAnalysisData = './dataExample/';
% pathToAnalysisData = '/Users/gomena/Research/EJBigData/EJ-2014-11-05-Processed/data005/';
pathToAnalysisData = '/Volumes/Analysis/2014-11-05-3/data005/';

neuronIds = [2418 2434 2913 2763]; % Enter the neurons that you want to 
                                   % analyze. *An elecResp file with each 
                                   % neuron number must exist*
patternNo = 1465;    % Stimulation pattern number. Should list just one.
recElecs  = [170 233 194 185]; % You can list any number. If you don't list
                               % any recElecs, then the script gets the recording
                               % electrodes from the elecResp file. If you
                               % list electrodes here, the script does not
                               % use the electrodes from elecResp. Instead,
                               % a function (prefElectrodes) is used to
                               % assign a main recording electrode to each
                               % of the neuronIds (chosen from this list of
                               % recElecs). The main recording electrode
                               % for each neuron is the electrode with the
                               % largest signal from the visual stim
                               % template. Other electrodes are used in the
                               % heuristics to improve the sigmoidal fit.
                               % Listing more recElecs is going to increase
                               % the algorithm runtime.


%Default minimum and maximum recording times. If the first part of the
%traces are not relevant, it could be worth trying make Tmin>1
Tmin     = 1;
Tmax     = 40;

% Load Templates from ElecResp. If recElecs are not supplied, then will be
% found from the 'goodElecs' field.
[templates, recElecs] = makeTemplatesFromElecResp(pathToAnalysisData,...
    patternNo,neuronIds,recElecs); %This function will output one electrode per neuron

% translate templates if by some reason there are weird template offsets.
templates = translateTemplate(templates,6,2,2); % The minimum of each 
                        % template should always align with sample point 10

%for each neuron, find the 'preFerred' electrode
prefElectrodes = PreferredElectrodes(templates); % Looks for the template with the maximum signal for each neuron
input.neuronInfo.prefElectrodes =   prefElectrodes;

%Parameters to find a Axonal activation authomatically. 
%if includeAxonBreakpoint=0 no axonal breakpoint will be included
input.params.load.findAxon.includeAxonBreakpoint  =   0; % For now, set to zero because the axonal breakpoint solution is not working
input.params.load.findAxon.numberOmit             =   4;
input.params.load.findAxon.lMovingAverage         =   1;
input.params.load.findAxon.TrangeAxon             =  [11 40];
%If one, look for decreases in the minimum of the eigenvalues of the
%variance (space) of the energy plot. If equal 2 looks for increases
input.params.load.findAxon.typeEigenvalue         =   2;
% cleanData = 1 erases the first trial of each movie (for cases where some
% of the trials did not result in electrical stimulation as intended)
input.params.load.cleanData                       =   1;
% if collapseTrialsSameCondition = 1 Collapse movies 2*j and 2*j-1
input.params.load.collapseTrialsSameCondition     =   0;

input.tracesInfo.Trange        =    [Tmin Tmax];
input.tracesInfo.recElecs      =    recElecs;
input.neuronInfo.neuronIds     =    neuronIds;
input.neuronInfo.templates     =    templates;

input = loadData(input,pathToAnalysisData,patternNo);


input = fillDefaultValues(input);

initial = Initialize(input);

[Gibbs GibbsNoDelete initial input Log] = doSpikeSortingElectricalArtifact(input,initial);
