function [Output]=DoSpikeSortingLargeScale(pathToPreparation,pathToEi,patternNo,neuronIds,params,varargin)
%Gonzalo Mena, 03/2016

Dif1=params.global.Dif1;
Dif2=params.global.Dif2;
positions=params.global.positions;
x=params.arrayInfo.x;
useStimElec = params.global.useStimElec;
useBundleAlg = params.bundle.useBundleAlg;

Tmax=params.global.Tmax;
nTrial=params.global.nTrial;

 

if(nargin==6)
   FoldersNames=varargin{1}; 
else
    
    
     
    dirs=dir(pathToPreparation);
    cont=1;
    
    for i=1:length(dirs)
    if(length(dirs(i).name)>=4)
        aux=find(dirs(i).name(1:4)=='data');
    if(length(aux)>0)
        
        FoldersNames{cont}=dirs(i).name;
        cont=cont+1;
    
    end
    end
    end
  
end

for f=1:length(FoldersNames)
    pathAux=[pathToPreparation FoldersNames{f}];
        dirs=dir(pathAux);
        
        
        
        for i=1:length(dirs)
        
        aux=find(dirs(i).name(1)=='p');
        if(length(aux)>0)
            
            if(isequal(dirs(i).name(2:end),num2str(patternNo)))
                path=pathAux;
            end
        end
        end
end

pathToAnalysisData=path;




if(params.global.sortData)
[TracesAll Art var0 listAmps listCurrents onset onsetC pval stimElecs]=loadTracesArtSort(pathToAnalysisData,patternNo,Tmax,nTrial,params);

else
    [TracesAll Art var0 listAmps listCurrents onset onsetC pval stimElecs]=loadTracesArt(pathToAnalysisData,patternNo,Tmax,nTrial,params);

end


if(length(stimElecs)>1)
    Output=1;
    disp('only one stimulating electrode supported')
    return
else
    stimElec=stimElecs;
end

[theta,rho] = cart2pol(positions(:,1)-positions(patternNo,1),positions(:,2)-positions(patternNo,2));


ind=setdiff([1:512],find(rho==0));

rho=rho(ind);

for j=1:length(listAmps)
for i=1:length(listAmps)
Dif3(i,j)=abs(listAmps(i)-listAmps(j));
end
end


Difs{2}=Dif2(ind,ind)/max(max(Dif2(ind,ind)));
Difs{3}=Dif3;
Difs{3}=Difs{3}/max(max(Difs{3}));
Dif1=Dif1/(max(max(Dif1)));
Difs{1}=Dif1;   

Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho/max(rho);
Diags{3}=listAmps'/max(listAmps);
     
factp=[0 3 6 9];

for k=1:3
[Ker KerD]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),[1 1 1]);
[a b]=eig(Ker);
Q{k}=a';
Qt{k}=a;
dL{k}=diag(b);
Kers{k}=Ker;
end



params.patternInfo.Q=Q;
params.patternInfo.Qt=Qt;
params.patternInfo.Kers=Kers;
params.patternInfo.dL=dL;
params.patternInfo.ind=ind;
params.patternInfo.Art =Art;
params.patternInfo.var0=var0;
params.patternInfo.rho=rho;
params.patternInfo.patternNo=patternNo;
params.patternInfo.listAmps=listAmps;
params.patternInfo.listCurrents=listCurrents;
params.patternInfo.Diags=Diags;
params.patternInfo.Difs=Difs;
params.patternInfo.stimElec=stimElec;


params.bundle.onsBundle=onsetC;
params.bundle.onset=onset;

templates=makeTemplatesFromEiShift(pathToEi, neuronIds,[1:512]);
params.neuronInfo.templates = templates;
if(~useStimElec&&(~useBundleAlg))
[spikes Log params]=SpikeSortingNoBundleNoStim(params,TracesAll);
elseif(~useStimElec&&(useBundleAlg))
[spikes Log params]=SpikeSortingBundleNoStim(params,TracesAll);    
elseif(useStimElec&&(~useBundleAlg))
    params = MakeStimKernels(params);
[spikes Log params]=SpikeSortingNoBundleStim(params,TracesAll);    
  elseif(useStimElec&&(useBundleAlg))
       params = MakeStimKernels(params);
[spikes Log params]=SpikeSortingBundleStim(params,TracesAll);    

end

if(useStimElec)
    Output.stimInfo.ActiveElectrodes=[1:512];
   Output.stimInfo.KersStim=params.patternInfo.KersStim;
    Output.stimInfo.breakpoints=params.patternInfo.breakpoints;
    Output.stimInfo.xStim=params.patternInfo.xStim;

else
    Output.stimInfo.ActiveElectrodes=ind;
end


Output.neuronInfo=params.neuronInfo;
Output.neuronInfo.neuronIds=neuronIds;
Output.neuronInfo.spikes=spikes;

Output.stimInfo.patternNo=patternNo;
Output.stimInfo.stimElec=stimElec;
Output.stimInfo.listAmps=listAmps;
Output.stimInfo.listCurrents=listCurrents;
Output.stimInfo.Art=params.patternInfo.Art;
Output.stimInfo.nTrials=params.patternInfo.nTrials;
Output.stimInfo.Kers=params.patternInfo.Kers;
Output.stimInfo.rho=rho;
Output.Log=Log;
Output.path.pathToEi=pathToEi;
Output.path.pathToPreparation=pathToPreparation;
Output.path.pathToAnalysisData=pathToAnalysisData;
Output.bundle=params.bundle;
Output.bundle.onsetC=onsetC;
Output.bundle.onset=onset;
Output.bundle.pvals=pval;
Output.stimInfo.useStimElec=useStimElec;
Output.stimInfo.Residuals=Res;
Output.arrayInfo=params.arrayInfo;
end

