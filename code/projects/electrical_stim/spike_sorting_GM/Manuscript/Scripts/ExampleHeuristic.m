clear input
    clear initial
    clear Gibbs
    clear GibbsNoDelete
    clear Log
    
    pathName = '/Users/gomena/Research/EJBigData/2015-04-14-0/data001/elecResp_n1699_p114.mat'
    temp = load(pathName);
    pathToAnalysisData =  '/Users/gomena/Research/EJBigData/2015-04-14-0/data001/'
    
    patternNo = temp.elecResp.stimInfo.patternNo;
    %pathToAnalysisData = temp.elecResp.names.data_path;
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
   % doSpikeSortingElectricalArtifact takes the input and initial structures as input and with that, executes the core
% part of the spike sorting algorithm: In the first place, it does Gibbs sampling iterations (iterating between
% spikes-latencies, Artifact, residual variances and logistic regression for spikes)
% until no changes are observed, or until a maximum number of iterations,
% specified in input.params.Gibbs.maxIter is exceeded (the choice of this number can be
% critical for the several neurons-electrodes context. Variables are stored in the Gibbs structure, and they suffer changes from
% iteration to iteration. After these first Gibbs sampler iterations, doSpikeSortingElectricalArtifact
% will apply sequentially the different heuristics (files that start with a H) in case some problem has been diagnosed
% (for example, lack of fit with the underlying logistic regression model)
% Almost all of these heuristics are based in re-estimating (extrapolating or interpolating)
% the artifact in conditions for which problems have been diagnosed, and doing more Gibbs iterations again, until
% no observed changes in spikes are observed (or until iterations exceed the maximum allowed).
% Nonetheless, one of the heuristics is more aggressive  and will try to delete ALL spikes in certain conditions, until this 
% delation leads to an increase in the residual variance. This strategy seems to work well in practice,
%and in the forst case, the solution withouth deletion will be part of the output as well.
% Notice that sequentiallness of Heuristics was done here only because it seemed that was the right order
%However, there are possibilities in exploring how the different Heuristics should succedd the others
%or even, if they could be included in a while loop that will be active until no more changes are needed.
%But for now, keeping things like this should provide reasonable results
%
%Input:        input and initial structure
%Output:       -same input and initial structure (this could be changed with no harm)
%              -Gibbs: a structure such that Gibbs.variables contains all relevant variables that are sampled 
%                      Also, Gibbs.params contains the paramters for the Gibbs sampler (as maximum number of iterations)
%                      And Gibbs.diagnostics contain information about logistic regression fit to spikes (not relevant now)
%              -GibbsNoDelete: same as Gibbs, but results reflect what happened before iterations
%              -Log: a structure such that Log(n) is the Log for neuron with index n, and
%               Log(n).Heuristic is a array of strings containing information about what changes were by the heuristics
%               (including interpolations/extrapolations) of artifact and deletion
%               Log(n).Deletion is a vector containing information about if it was or no deletion for neuron n at the different breakpoints ranges.
%Important:     Essentially, spike sorting solutions correspond to Gibbs.variables.spikes and Gibbs.variables.latencies (the same with GibbsNoDelete, in case
%               the experimenter thinks they can give useful information



%%
%load relevant variables and create Gibbs structure
data           = input.tracesInfo.data;
dataVecJ       = input.tracesInfo.dataVecJ;
templates      = input.neuronInfo.templates;
TfindRangeRel  = input.params.initial.TfindRangeRel;
TfindRel       = [TfindRangeRel(1):TfindRangeRel(2)];
TfindRel0      = [0 TfindRel]; 
neuronIds      = input.neuronInfo.neuronIds;

E = input.tracesInfo.E;
T = input.tracesInfo.T;
I = input.tracesInfo.I;
J = input.tracesInfo.J;

nNeurons = input.neuronInfo.nNeurons;

for n=1:nNeurons
    Log(n).params.contLogHeuristic = 0;
    for j=1:J
        Gibbs.variables.spikes{n}(j,1:I(j))=100000;
    end
end

lengthDataVec   = sum(I)*E*T;



K = makeToeplitz(templates,TfindRel0,T);

for n=1:nNeurons
    
  
    indn=zeros(1,nNeurons);
    indn(n)=1;
    
    for t=1:length(TfindRel0)

        indt       = zeros(1,length(TfindRel0));
        indt(t)    = 1;
        indnt      = kron(indn,indt);
        Kn{n}(:,t) = K*indnt';
        end
end


Gibbs.params.I = I;
Gibbs.params.E = E;
Gibbs.params.J = J;
Gibbs.params.T = T;
Gibbs.params.nNeurons = nNeurons;
Gibbs.params.maxIterGibbs = input.params.maxIterGibbs;

Gibbs.params.Kn                     = Kn;
Gibbs.params.X                      = initial.params.X;
Gibbs.params.Xj                     = initial.params.Xj;

for e=1:E
    Gibbs.params.matricesReg(e).Prods = initial.params.matricesReg(e).Prods;
end

Gibbs.params.lambdas                = initial.params.lambdas;
Gibbs.params.Lambdas                = initial.params.Lambdas;
Gibbs.params.LambdasInv             = initial.params.LambdasInv;
Gibbs.params.TfindRel0              = TfindRel0;
Gibbs.params.Tfind0                 = [0 input.params.initial.TfindRange(1):input.params.initial.TfindRange(2)];
Gibbs.params.dataVecJ               = dataVecJ;
Gibbs.params.a0                     = initial.params.a0;
Gibbs.params.b0                     = initial.params.b0;
Gibbs.params.lambdaLogReg           = initial.params.lambdaLogReg;
Gibbs.params.alphaLogReg            = initial.params.alphaLogReg;
Gibbs.params.Tdivide                = input.params.initial.Tdivide;
Gibbs.params.thresLogistic          = input.params.Gibbs.thresLogistic;

Gibbs.variables.sigma            = initial.sigma;
Gibbs.variables.Artifact         = initial.Artifact;
Gibbs.variables.Probs            = initial.Probs;
Gibbs.variables.ActionPotentials = initial.ActionPotentials;  
Gibbs.variables.Beta             = zeros(size(Gibbs.params.Lambdas{1},1),E);     
%%

%Do first Gibbs sampler iterations
Gibbs = GibbsSamplerSpikesArtifact(Gibbs);
S{1}=Gibbs
[Gibbs Log] = HdeleteSpikesBeginning(input,Gibbs,Log,1);
S{2}=Gibbs;


cont=2;
type        = input.params.Heuristic.ResampleAlltimes;
E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
nNeurons    = Gibbs.params.nNeurons;

thres          = input.params.Heuristic.LogisticRegPValueThres;
maxCondsChange = input.params.Heuristic.maxCondsChange;



flags = zeros(nNeurons,1);

for n=1:nNeurons
    contLog(n) = Log(n).params.contLogHeuristic;
    CondsChangeold{n}=[];
    nCondsChange(n)=1;
end

while(prod(flags)==0)

    p2          = Gibbs.diagnostics.Logisticsp2GoF;
    condsSorted = Gibbs.diagnostics.condSortError;
        
    for n=1:nNeurons

        if(p2(n)>thres)
            
            CondsChange{n}=[];
            CondsChangeold{n}=CondsChange{n};
            flags(n)=1;
            continue
    
        end


        CondsChange{n}=condsSorted(n,1:min(nCondsChange(n),maxCondsChange));

        if(length(CondsChange{n})==maxCondsChange)
            random          = unidrnd(maxCondsChange-1);
            nCondsChange(n) = random;
            randSample      = sort(randsample(maxCondsChange,random))';
            CondsChange{n}  = condsSorted(n,randSample);
            flags(n)=1;
    
            continue
        elseif(isequal(unique(CondsChange{n}),unique(CondsChangeold{n})))
            nCondsChange(n)=nCondsChange(n)+1;
            CondsChange{n}=condsSorted(n,1:min(nCondsChange(n),maxCondsChange));

        else
            nCondsChange(n)=1;
            CondsChange{n}=condsSorted(n,1:min(nCondsChange(n),maxCondsChange));

        end


    CondsChangeold{n}=CondsChange{n};

    end

    CondsChange    = CondsNeuron2Electrode(input,CondsChange);
    
    for n=1:nNeurons
        if(~isempty(CondsChange{n}))
            contLog(n)=contLog(n)+1;
            Log(n).Heuristic{contLog(n)}=['Logistic regression fit is poor at Conditions ' num2str(CondsChange{n}) ' p = ' num2str(Gibbs.diagnostics.Logisticsp2GoF(n))];
            Log(n).params.contLogHeuristic = contLog(n);
        end
    end
    ArtifactNew=[];

    for e = 1:E
    
        if(type ==1)
            
            Gibbs       = SampleConditionalArtifact(Gibbs,e,CondsChange{e},setdiff([1:J],CondsChange{e}));
        else
            Gibbs       = SampleConditionalArtifactT(Gibbs,e,CondsChange{e},setdiff([1:J],CondsChange{e})); 
        end
        ArtifactNew = [ArtifactNew Gibbs.variables.ArtifactE{e}];
   
   
    end
    Artifact                 = ArtifactNew;
    Gibbs.variables.Artifact = ArtifactNew;
    Gibbs                    = GibbsSamplerSpikesArtifact(Gibbs);
    cont=cont+1;
    S{cont}=Gibbs;

end


cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/Scripts')
amps=abs(input.stimInfo.listAmps)';

%%

cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/FiguresMethods/AlgExample')
colors=colormap(jet);
time=[1:40]/20;
for m=1:size(initial.Artifact,1)
    
    plot(time,initial.Artifact(m,:),'linewidth',2,'color',colors(floor(m*64/39),:))
    hold on
    xlabel('time (ms)','fontsize',20)
    ylabel('recorded daqs','fontsize',20)
    set(gca,'Fontsize',14)
    title('Initial Artifact','fontsize',20)
    
end
print('A0','-djpeg')
    print('A0','-depsc2')
listCond=[length(S)-3:length(S)-1];


condBlack{1}=[22];
condBlack{2}=[23];
condBlack{3}=[];
legends={'Artifact after Gibbs sampling','Artifact after Gibbs sampling','Artifact after heuristic'};
for c=1:length(listCond)
    c2=listCond(c);
    figure
    for m=1:size(initial.Artifact,1)
        
        plot(time,S{c2}.variables.ArtifactE{1}(m,:),'linewidth',2,'color',colors(floor(m*64/39),:))
        hold on
        xlabel('time (ms)','fontsize',20)
        ylabel('recorded daqs','fontsize',20)
        set(gca,'Fontsize',14)
        title(legends{c},'fontsize',20)
        
    end
    for b=1:length(condBlack{c})
         plot(time,S{c2}.variables.ArtifactE{1}(condBlack{c}(b),:),'--','linewidth',5,'color',colors(floor(condBlack{c}(b)*64/39),:))
    end
    print(['A' num2str(c)'],'-djpeg')
        print(['A' num2str(c)'],'-depsc2')
end



for ii = 1:size(elecResp.analysis.latencies,1)
    humanLats = elecResp.analysis.latencies{ii}; 
    humanLats = humanLats(humanLats>0); 
    humanLatencies(ii,1) = prctile(humanLats,25);
    humanLatencies(ii,2) = prctile(humanLats,50);
    humanLatencies(ii,3) = prctile(humanLats,75);
end

humanProb=elecResp.analysis.successRates';


legends={'Activation curves after Gibbs sampling','Activation curves after Gibbs sampling','Activation curves after after heuristic'};

for c=1:length(listCond)
    figure
    spikeProbs=nansum(S{listCond(c)}.variables.spikes{1}')./input.tracesInfo.I;
    spikeLogProbs=S{listCond(c)}.variables.ProbRobust;
    
    
    hum=plot(amps,humanProb,'ro-','LineWidth',1.5);
    hold on; alg = plot(amps,spikeProbs,'^-','linewidth',1.5);
    logRegAlg = plot(amps,spikeLogProbs,'linewidth',1.5,'color','black');
    grid('on')
    xlabel('Stimulus amplitude (\muA)','fontsize',20)
    ylabel('Spike Probability','fontsize',20)
    set(gca,'Fontsize',14)
    title(legends{c},'fontsize',20)
     for b=1:length(condBlack{c})
         scatter(amps(condBlack{c}(b)),1,100,'green','filled')
    end
    %legend([hum alg logRegAlg],'Human','Algorithm','Logistic Regression')
        print(['Act' num2str(c)'],'-djpeg')
        print(['Act' num2str(c)'],'-depsc2')
end


legends={'Activation curves after Gibbs sampling','Activation curves after Gibbs sampling','Activation curves after after heuristic'};


    figure
    spikeProbs=initial.Probs;
    %spikeLogProbs=S{listCond(c)}.variables.ProbRobust;
    
    
    hum=plot(amps,humanProb,'ro-','LineWidth',1.5);
    hold on; alg = plot(amps,spikeProbs,'^-','linewidth',1.5);
    %logRegAlg = plot(amps,spikeLogProbs,'linewidth',1.5,'color','black');
    grid('on')
    xlabel('Stimulus amplitude (\muA)','fontsize',20)
    ylabel('Spike Probability','fontsize',20)
    set(gca,'Fontsize',14)
    title('Initial activation curves','fontsize',20)
    %legend([hum alg logRegAlg],'Human','Algorithm','Logistic Regression')
        print(['Act0'],'-djpeg')
        print(['Act0'],'-depsc2')




