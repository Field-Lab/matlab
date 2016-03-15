function [spikes Log params]=SpikeSortingBundleNoStim(params,TracesAll)


Kers=params.patternInfo.Kers;
Q=params.patternInfo.Q;
Qt=params.patternInfo.Qt;
dL=params.patternInfo.dL;

ind=params.patternInfo.ind;
templates=params.neuronInfo.templates;
Art=params.patternInfo.Art;
Difs=params.patternInfo.Difs;
Diags=params.patternInfo.Diags;
var0=params.patternInfo.var0;
patternNo=params.patternInfo.patternNo;

thresEI=params.global.thresEI;
Tmax=params.global.Tmax;
tarray=params.global.tarray;
maxIter=params.global.maxIter*length(templates);
options= params.global.options;

cutBundle=params.bundle.cutBundle;
nVec=params.bundle.nVec;
updateFreq=params.bundle.updateFreq;



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
    els=union(b(ind2),els);
end


for n=1:length(templates)
    for t=1:length(tarray)
        [ActionPotential]=makeActionPotential(n,tarray(t),templates,Tmax);
        
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
params.patternInfo.nTrials(i)=trialI;
cont=1;
while(flag==1&&cont<=maxIter)
    
    clear times
    
    ArtF=FilterArtifactLocal(Kers,Art(1,ind,1:Tmax),[x log(var0)],i,ind,Q,Qt,krondiaginv);
    
    
    AA0=reshape(ArtF(i,els,:),Tmax*length(ind(els)),1);
    
    r=randsample(length(templates),length(templates));
    
    
    
    TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
    
    for n=1:length(templates)
        
        AA=reshape(TracesResidual(1:trialI,ind(els),1:Tmax),trialI,Tmax*length(ind(els)))'-repmat(AA0,1,trialI);
        
        
        corrs=-2*AA'*Kn{r(n)}+repmat(nansum(Kn{r(n)}.^2),trialI,1);
        [mins tmax]=min(corrs');
        times(r(n),:)=tmax;
        
        TracesResidual(:,ind,:)=TracesResidual(:,ind,:)-Knn{r(n)}(tmax,:,:);
        
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
        spikes{n}(1,1:trialI)=tarray(times(n,:));
    end
    cont=cont+1;
end
Log(1)=cont;


xold=x;
for i=2:maxCond
    
    
    if(i>=params.bundle.onsBundle)
        
        if(i==params.bundle.onsBundle)
            [Res]=ResidualsElectrodeSimple(Art(params.bundle.onsBundle:end,:,:),patternNo,[1:Tmax]);
            [a bb c]=svd(Res(:,ind));
            
            v1=max(a(:,1:nVec)*bb(1:nVec,1:nVec)*c(:,1:nVec)',0);
        end
        
        
        Difs{3}=0;
        Diags{3}=1;
        
        DiagsBundle=Diags;
        DiagsBundle{2}=v1(i-params.bundle.onsBundle+1,:)';
        if(mod(i-params.bundle.onsBundle,updateFreq)==0)
        f11=@(Art,x)logDetKron(Art(i,ind,:),[xold(1:9) x log(var0)],Difs,10,[1 4 1],DiagsBundle,[3 3 3]);
        
        g11=@(x)f11(Art,x);
        x0=0;
        x11 = fminunc(g11,x0,options);
        x(end)=x11;
        end
        
        k=2;
        [Ker KerD]=evalKernels(Difs{k},DiagsBundle{k},[xold(4:6) ],4);
        
         KersNew{k}=Ker;
        
        
        [a b]=eig(KersNew{k});
        Q{k}=a';
        Qt{k}=a;
        dL{k}=diag(b);
        Kers{k}=KersNew{k};
        
        
        krondiag0=kron(dL{1},dL{2});
        
        
        
        end
    
        
    krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
    
    trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
    params.patternInfo.nTrials(i)=trialI;
    [Apred]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,i,ArtF,x,var0);
    
    flag=1;
    cont=1;
    while(flag==1&&cont<=maxIter)
        
        
        clear times
        
        AA0=reshape(Apred(:,els,:),Tmax*length(ind(els)),1);
        
        
        r=randsample(length(templates),length(templates));
        
        
        
        TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
        for n=1:length(templates)
            
            AA=reshape(TracesResidual(1:trialI,ind(els),1:Tmax),trialI,Tmax*length(ind(els)))'-repmat(AA0,1,trialI);
            
            
            corrs=-2*AA'*Kn{r(n)}+repmat(nansum(Kn{r(n)}.^2),trialI,1);
            [mins tmax]=min(corrs');
            times(r(n),:)=tmax;
            TracesResidual(:,ind,:)=TracesResidual(:,ind,:)-Knn{r(n)}(tmax,:,:);
            
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
            spikes{n}(i,1:trialI)=tarray(times(n,:));
        end
        cont=cont+1;
    end
    
    Log(i)=cont;
end
params.patternInfo.Art=Art(1:maxCond,:,:);