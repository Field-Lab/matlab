function [spikes Log params]=SpikeSortingNoBundleStim(params,TracesAll)
%Gonzalo Mena, 3/2016

Kers=params.patternInfo.Kers;
KersStim=params.patternInfo.KersStim;
dLStim=params.patternInfo.dLStim;
QStim=params.patternInfo.QStim;;
QtStim=params.patternInfo.QtStim;
Q=params.patternInfo.Q;
Qt=params.patternInfo.Qt;
dL=params.patternInfo.dL;
ind=params.patternInfo.ind;
templates=params.neuronInfo.templates;
Art=params.patternInfo.Art;
var0=params.patternInfo.var0;
patternNo=params.patternInfo.patternNo;
listCurrents=params.patternInfo.listCurrents;

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
xStim=params.patternInfo.xStim;


els=[];
for n=1:length(templates)
    templates{n}=templates{n}(:,:);
    spikes{n}=NaN*zeros(maxCond,size(TracesAll,2));
    [a b]=sort(max(abs(templates{n}')),'descend');
    ind2=find(a>thresEI);
    params.neuronInfo.ActiveElectrodes{n}=b(ind2);
    els=union(b(ind2),els);
end

params.neuronInfo.ActiveElectrodesAll=els;


elExtra=zeros(size(Art,1),1);
elExtra(logical([0;diff(listCurrents)>0]))=patternNo;

indpat=find(max(elExtra)==els);

for i=1:length(elExtra);
    if(elExtra(i)==patternNo)
         indels{i}=setdiff(1:Tmax*length(els),indpat:length(els):Tmax*length(els));
    else
       
    indels{i}=1:Tmax*length(els);
    end
end



for n=1:length(templates)
    for t=1:length(tarray)
        [ActionPotential]=makeActionPotential(n,tarray(t),templates,Tmax);
        
        Knn{n}(t,:,:)=ActionPotential(:,:);
        Kn{n}(:,t)=reshape(ActionPotential(els,:),Tmax*length(els),1)';
    end
end


flag=1;

krondiag0=1;
for k=1:2
    krondiag0=kron(krondiag0,dL{k});
end


krondiag0Stim=1;

    krondiag0Stim=kron(krondiag0Stim,dLStim{1});

i=1;

krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
krondiaginvStim=(exp(xStim(end))*krondiag0Stim*KersStim{2}(i,i)+var0).^(-1);


els2=setdiff(els,elExtra(i));

trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
params.patternInfo.nTrials(i)=trialI;
cont=1;
while(flag==1&&cont<=maxIter)
    
    clear times
    
    ArtF(:,ind,:)=FilterArtifactLocal(Kers,Art(1,ind,1:Tmax),[x log(var0)],i,ind,Q,Qt,krondiaginv);
    ArtF(:,patternNo,:)=FilterArtifactLocalStim(KersStim,squeeze(Art(1,patternNo,1:Tmax))',[xStim log(var0)],i,QStim,QtStim,krondiaginvStim);

    
    AA0=reshape(ArtF(i,els2,:),Tmax*length(els2),1);
    
    r=randsample(length(templates),length(templates));
    
    
    
    TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
    
    for n=1:length(templates)
        
        AA=reshape(TracesResidual(1:trialI,els,1:Tmax),trialI,Tmax*length(els))'-repmat(AA0,1,trialI);
        
        
        corrs=-2*AA'*Kn{r(n)}(indels{i},:)+repmat(nansum(Kn{r(n)}(indels{i},:).^2),trialI,1);
        [mins tmax]=min(corrs');
        times(r(n),:)=tmax;
        
        TracesResidual(:,:,:)=TracesResidual(:,:,:)-Knn{r(n)}(tmax,:,:);
        
    end
    
    Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
    
ArtF(1,ind,:)=FilterArtifactLocal(Kers,Art(1,ind,1:Tmax),[x log(var0)],1,ind,Q,Qt,krondiaginv);
ArtF(1,patternNo,:)=FilterArtifactLocalStim(KersStim,squeeze(Art(1,patternNo,1:Tmax))',[xStim log(var0)],1,QStim,QtStim,krondiaginvStim);
  
    flag2=ones(length(templates),1);
    for n=1:length(templates)
        
        if(nansum(tarray(times(n,:))==spikes{n}(1,1:trialI))==trialI)
            flag2(n)=0;
        end
    end
    flag=max(flag2);
    for n=1:length(templates)
        spikes{n}(1,1:trialI)=tarray(times(n,:));
        
    end
    cont=cont+1;
end
Log.Iter(1)=cont;



for i=2:maxCond
    
    
    
    
    krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
    krondiaginvStim=(exp(xStim(end))*krondiag0Stim*KersStim{2}(i,i)+var0).^(-1);

    trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
    params.patternInfo.nTrials(i)=trialI;
    
[Apred1]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,i,ArtF(:,ind,:),x,var0);   
[Apred2]=ExtrapolateArtifactCondStim(KersStim,QStim,QtStim,dLStim,i,squeeze(ArtF(:,patternNo,:)),xStim,var0);   
Apred(:,ind,:)=Apred1;
Apred(:,patternNo,:)=Apred2;

    flag=1;
    cont=1;
    while(flag==1&&cont<=maxIter)
        
        
        clear times
        
        AA0=reshape(Apred(:,els,:),Tmax*length(ind(els)),1);
        
        
        r=randsample(length(templates),length(templates));
        
        
        
        TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
        for n=1:length(templates)
            
            AA=reshape(TracesResidual(1:trialI,ind(els),1:Tmax),trialI,Tmax*length(ind(els)))'-repmat(AA0,1,trialI);
            
            
            corrs=-2*AA'*Kn{r(n)}(indels{i},:)+repmat(nansum(Kn{r(n)}(indels{i},:).^2),trialI,1);
            [mins tmax]=min(corrs');
            times(r(n),:)=tmax;
            TracesResidual(:,:,:)=TracesResidual(:,:,:)-Knn{r(n)}(tmax,:,:);
            
        end
        
        
        Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
        
        ArtF(i,ind,:)=FilterArtifactLocal(Kers,Art(1:i,ind,:),[x log(var0)],i,ind,Q,Qt,krondiaginv);
        ArtF(i,patternNo,:)=FilterArtifactLocalStim(KersStim,squeeze(Art(1:i,patternNo,:)),[xStim log(var0)],i,QStim,QtStim,krondiaginvStim);

        
        Apred=ArtF(i,:,:);
        
        
        flag2=ones(length(templates),1);
        for n=1:length(templates)
            
            if(nansum(tarray(times(n,:))==spikes{n}(i,1:trialI))==trialI)
                flag2(n)=0;
            end
        end
        flag=max(flag2);
        for n=1:length(templates)
            spikes{n}(i,1:trialI)=tarray(times(n,:));
            
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
