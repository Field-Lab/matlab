function [spikes Log params]=SpikeSortingNoBundleNoStim(params,TracesAll)
%Gonzalo Mena, 3/2016

Kers=params.patternInfo.Kers;
Q=params.patternInfo.Q;
Qt=params.patternInfo.Qt;
dL=params.patternInfo.dL;
ind=params.patternInfo.ind;
templates=params.neuronInfo.templates;
Art=params.patternInfo.Art;
var0=params.patternInfo.var0;

thresEI=params.global.thresEI;
Tmax=params.global.Tmax;
tarray=params.global.tarray;
maxIter=params.global.maxIter*length(templates);
cutBundle=params.bundle.cutBundle;
contMessage=1;

if(cutBundle==1)
    maxCond=params.bundle.onsBundle-1;
else
    maxCond=size(TracesAll,1);
end


x=params.arrayInfo.x;

els=[];
for n=1:length(templates)
    templates{n}=templates{n}(ind,:);
    spikes{n}=NaN*zeros(maxCond,size(TracesAll,2));
    [a b]=sort(max(abs(templates{n}')),'descend');
    ind2=find(a>thresEI);
    if(isempty(ind2))
        ind2=1;
    else
    
    end
    
    params.neuronInfo.ActiveElectrodes{n}=ind(b(ind2));
    els=union(b(ind2),els);
end

tarray2=setdiff(tarray,0);
params.neuronInfo.ActiveElectrodesAll=ind(els);
KnAll=zeros(length(els)*Tmax,1);
indNeurons=[0 kron([1:length(templates)],ones(1,length(tarray2)))];
indTimes=[0 kron(ones(1,length(templates)),tarray2)];
indTimesArray=[0 kron(ones(1,length(templates)),[1:length(tarray2)])];

for n=1:length(templates)
    for t=1:length(tarray2)
        [ActionPotential]=makeActionPotential(n,tarray2(t),templates,Tmax);
        
        Knn(n,t,:,:)=ActionPotential(:,:);
        Aux(:,t)=reshape(ActionPotential(els,:),Tmax*length(ind(els)),1)';
       
    end
     KnAll=[KnAll Aux];
end

KnnReshaped=reshape(Knn,size(Knn,1)*size(Knn,2),size(Knn,3),size(Knn,4));
        
flag=1;

krondiag0=1;
for k=1:2
    krondiag0=kron(krondiag0,dL{k});
end
i=1;

krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);


trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
params.patternInfo.nTrials(i)=trialI;
cont=1;
    
while(flag==1&&cont<=maxIter)
    
    clear times
    
    ArtF=FilterArtifactLocal(Kers,Art(1,ind,1:Tmax),[x log(var0)],i,ind,Q,Qt,krondiaginv);
    
    
    AA0=reshape(ArtF(i,els,:),Tmax*length(ind(els)),1);
    
   
    TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
    
    indTrial=[1:trialI];
    contFindNeurons=1;
    times=zeros(length(templates),trialI);
    while(length(indTrial)>0)
        AA=reshape(TracesResidual(indTrial,ind(els),1:Tmax),length(indTrial),Tmax*length(ind(els)))'-repmat(AA0,1,length(indTrial));
        
        
        corrs=-2*AA'*KnAll+repmat(nansum(KnAll.^2),length(indTrial),1);
        if(contFindNeurons>1)
        for trial=1:length(indTrial)
        corrs(trial,find(indNeurons==indNeurons(tmax(indTrialRel(trial)))))=NaN;
        
        end
        end
        [mins tmax]=nanmin(corrs');
        
        indTrialRel=find(tmax>1);
        indTrial=indTrial(indTrialRel);
         idx = sub2ind(size(times), indNeurons(tmax(indTrialRel)),indTrialRel);
        times(idx)=indTimes(tmax(indTrialRel));
   
        idxSubtract = sub2ind(size(squeeze(Knn(:,:,1,1))), indNeurons(tmax(indTrialRel)),indTimesArray(tmax(indTrialRel)));
    
        TracesResidual(indTrial,ind,:)=TracesResidual(indTrial,ind,:)-KnnReshaped(idxSubtract,:,:);
       contFindNeurons=contFindNeurons+1;
    end
    
    Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
    
    ArtF=FilterArtifactLocal(Kers,Art(1,ind,1:Tmax),[x log(var0)],1,ind,Q,Qt,krondiaginv);
    
    flag2=ones(length(templates),1);
    for n=1:length(templates)
        
        if(nansum(times(n,:)==spikes{n}(1,1:trialI))==trialI)
            flag2(n)=0;
        end
    end
    flag=max(flag2);
    for n=1:length(templates)
        spikes{n}(1,1:trialI)=times(n,:);
        
    end
    cont=cont+1;
end
Log.Iter(1)=cont;



for i=2:maxCond
    
    
i  

    
    krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
    
    trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
    params.patternInfo.nTrials(i)=trialI;
    
    [Apred]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,i,ArtF,x,var0);
    
    flag=1;
    cont=1;
    while(flag==1&&cont<=maxIter)
        
        
        clear times
        
        AA0=reshape(Apred(:,els,:),Tmax*length(ind(els)),1);
        
        
  
    TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
    
    indTrial=[1:trialI];
    contFindNeurons=1;
    times=zeros(length(templates),trialI);
    while(length(indTrial)>0)
        AA=reshape(TracesResidual(indTrial,ind(els),1:Tmax),length(indTrial),Tmax*length(ind(els)))'-repmat(AA0,1,length(indTrial));
        
        
        corrs=-2*AA'*KnAll+repmat(nansum(KnAll.^2),length(indTrial),1);
        if(contFindNeurons>1)
        for trial=1:length(indTrial)
        corrs(trial,find(indNeurons==indNeurons(tmax(indTrialRel(trial)))))=NaN;
        
        end
        end
        [mins tmax]=nanmin(corrs');
        
        indTrialRel=find(tmax>1);
        indTrial=indTrial(indTrialRel);
         idx = sub2ind(size(times), indNeurons(tmax(indTrialRel)),indTrialRel);
        times(idx)=indTimes(tmax(indTrialRel));
   
        idxSubtract = sub2ind(size(squeeze(Knn(:,:,1,1))), indNeurons(tmax(indTrialRel)),indTimesArray(tmax(indTrialRel)));
    
        TracesResidual(indTrial,ind,:)=TracesResidual(indTrial,ind,:)-KnnReshaped(idxSubtract,:,:);
       contFindNeurons=contFindNeurons+1;
    end
   
    
    
        Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
        
        ArtF(i,:,:)=FilterArtifactLocal(Kers,Art(1:i,ind,:),[x log(var0)],i,ind,Q,Qt,krondiaginv);
        
        
        Apred=ArtF(i,:,:);
        
        
        flag2=ones(length(templates),1);
        for n=1:length(templates)
            
            if(nansum(times(n,:)==spikes{n}(i,1:trialI))==trialI)
                flag2(n)=0;
            end
        end
        flag=max(flag2);
        for n=1:length(templates)
            spikes{n}(i,1:trialI)=times(n,:);
            
        end
        cont=cont+1;
        if(cont==maxIter)
             Log.Message{contMessage}=['Maximum number of iterations exceeded at conditon ' num2str(i) ];
        contMessage=contMessage+1;
        end
    end
    
    Log.Iter(i)=cont-1;
end
params.patternInfo.Art=Art(1:maxCond,:,:);
