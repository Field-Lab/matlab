% Example bug




clear input
clear initial
clear Gibbs
clear GibbsNoDelete
clear Log

temp = load('/Users/gomena/Research/EJBigData/Datasetsvisitjun15/Debugging/elecResp_n2796_p138.mat');
patternNo = temp.elecResp.stimInfo.patternNo;
pathToAnalysisData = '/Users/gomena/Research/EJBigData/Datasetsvisitjun15/Debugging/';
neuronIds = temp.elecResp.cells.main;


% analyze. *An elecResp file with each
% neuron number must exist*
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
input.tracesInfo.Trange        =    [Tmin Tmax];


% Load Templates from ElecResp. If recElecs are not supplied, then will be
% found from the 'goodElecs' field.
[templates, recElecs] = makeTemplatesFromElecResp(pathToAnalysisData,...
    patternNo,neuronIds,1); %This function will output one electrode per neuron
if(size(templates{1},2)==81) %Usually, when templates have length 81, they are aligned
    %To spike onset at time ~20. The following
    %line re-aligns them to onset at time 10
    [templates, recElecs] = makeTemplatesFromElecResp(pathToAnalysisData,...
        patternNo,neuronIds,11);
end

%translate templates if by some reason there are weird template offsets.
%templates = translateTemplate(templates,Translate(p),1,1); % The minimum of each
% template should always align with sample point 10

%for each neuron, find the 'preFerred' electrode
prefElectrodes = PreferredElectrodes(templates); % Looks for the template with the maximum signal for each neuron
input.neuronInfo.prefElectrodes =   prefElectrodes; %Define preferred electrodes
input.tracesInfo.recElecs      =    recElecs; %define recording electrodes
input.neuronInfo.neuronIds     =    neuronIds; %define neuron Ids
input.neuronInfo.templates     =    templates; %define templates

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
input.params.load.collapseTrialsSameCondition     =   1;


%loads data from movie files
input = loadData(input,pathToAnalysisData,patternNo);

%Fill default values for input
input = fillDefaultValues(input);

%any change in the default values should be done at this point;
%input.params.initial.degPolRule{e}=[1 2 3 50]; example of a redefinition of a default value

%finds initial values for the artifact and rest of variables (via convex
%relaxation with cvx-Mosek)
initial = Initialize(input);

%The Core of the algorithm: spike sorting via gibbs sampler + heuristics.
%solutions are stored in Gibbs Structure
[Gibbs GibbsNoDelete initial input Log] = doSpikeSortingElectricalArtifact(input,initial);

Output_example = OutputResults(input,Gibbs,Log);
