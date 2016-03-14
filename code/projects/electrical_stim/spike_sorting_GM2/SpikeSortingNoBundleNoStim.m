function [spikes Log]=SpikeSortingNoBundleNoStim(params)


Kers=params.patternInfo.Kers;
Q=params.patternInfo.Q;
Qt=params.patternInfo.Qt;
dL=params.patternInfo.dL;
TracesAll=params.patternInfo.TracesAll;
ind=params.patternInfo.ind;
template=params.neuronInfo.template;
Art=params.patternInfo.Art0;


thresEI=params.global.thresEI;
Tmax=params.global.Tmax;
tarray=params.global.tarray;
maxIter=params.global.maxIter;
cutBundle=params.bundle.cutBundle;
nvect=params.bundle.nVec;

if(cutBundle==1)
    maxCond=params.bundle.onsCond;
else
    maxCond=size(TracesAll,1);
end


x=params.arrayInfo.x;



els=[];
for n=1:length(template)
    spikes{n}=NaN*zeros(size(TracesAll,1),size(TracesAll,2));
    [a b]=sort(max(abs(template{n}')),'descend');
    ind2=find(a>thresEI);
    els=union(b(ind2),els);
end


for n=1:length(template)
    for t=1:length(tarray)
        [ActionPotential]=makeActionPotential(n,tarray(t),template,Tmax);
        
        Knn{n}(t,:,:)=ActionPotential(:,:);
        Kn{n}(:,t)=reshape(ActionPotential(els,:),Tmax*length(ind(els)),1)';
    end
end


flag=1;

krondiag0=1;
for k=1:2
    krondiag0=kron(krondiag0,dL{k});
end
i=1;

krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);


trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
cont=1;
while(flag==1&&cont<=maxIter)
    
    clear times
    
    ArtF=FilterArtifactLocal(Kers,Art(1,ind,1:Tmax),[x log(var0)],i,ind,Q,Qt,krondiaginv);
    
    
    AA0=reshape(ArtF(i,els,:),Tmax*length(ind(els)),1);
    
    r=randsample(length(template),length(template));
    
    
    
    TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
    
    for n=1:length(template)
        
        AA=reshape(TracesResidual(1:trialI,ind(els),1:Tmax),trialI,Tmax*length(ind(els)))'-repmat(AA0,1,trialI);
        
        
        corrs=-2*AA'*Kn{r(n)}+repmat(nansum(Kn{r(n)}.^2),trialI,1);
        [mins tmax]=min(corrs');
        times(r(n),:)=tmax;
        
        TracesResidual(:,ind,:)=TracesResidual(:,ind,:)-Knn{r(n)}(tmax,:,:);
        
    end
    
    Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
    
    ArtF=FilterArtifactLocal(Kers,Art(1,ind,1:Tmax),[x log(var0)],1,ind,Q,Qt,krondiaginv);
    
    flag2=ones(length(template),1);
    for n=1:length(template)
        
        if(nansum(times(n,:)==spikes{n}(1,1:trialI))==trialI)
            flag2(n)=0;
        end
    end
    flag=max(flag2);
    for n=1:length(template)
        spikes{n}(1,1:trialI)=tarray(times(n,:));
    end
    cont=cont+1;
end
Log(1)=cont;



for i=2:maxCond
    
    
    
    
    krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
    
    trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
    [Apred]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,i,ArtF,x,var0);
    
    flag=1;
    cont=1;
    while(flag==1&&cont<=maxIter)
        
        
        clear times
        
        AA0=reshape(Apred(:,els,:),Tmax*length(ind(els)),1);
        
        
        r=randsample(length(template),length(template));
        
        
        
        TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
        for n=1:length(template)
            
            AA=reshape(TracesResidual(1:trialI,ind(els),1:Tmax),trialI,Tmax*length(ind(els)))'-repmat(AA0,1,trialI);
            
            
            corrs=-2*AA'*Kn{r(n)}+repmat(nansum(Kn{r(n)}.^2),trialI,1);
            [mins tmax]=min(corrs');
            times(r(n),:)=tmax;
            TracesResidual(:,ind,:)=TracesResidual(:,ind,:)-Knn{r(n)}(tmax,:,:);
            
        end
        
        
        Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
        
        ArtF(i,:,:)=FilterArtifactLocal(Kers,Art(1:i,ind,:),[x log(var0)],i,ind,Q,Qt,krondiaginv);
        
        
        Apred=ArtF(i,:,:);
        
        
        flag2=ones(length(template),1);
        for n=1:length(template)
            
            if(nansum(times(n,:)==spikes{n}(i,1:trialI))==trialI)
                flag2(n)=0;
            end
        end
        flag=max(flag2);
        for n=1:length(template)
            spikes{n}(i,1:trialI)=tarray(times(n,:));
        end
        cont=cont+1;
    end
    
    Log(i)=cont;
end
