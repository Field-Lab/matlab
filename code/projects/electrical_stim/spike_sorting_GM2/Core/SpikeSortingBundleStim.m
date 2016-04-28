function [spikes Log params]=SpikeSortingNoBundleStim(params,TracesAll)
%Gonzalo Mena, 3/2016

Kers=params.patternInfo.Kers;
KersSti=params.patternInfo.KersStim;
dLSti=params.patternInfo.dLStim;
QSti=params.patternInfo.QStim;
QtSti=params.patternInfo.QtStim;
options=params.global.options;
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

stimElec=params.patternInfo.stimElec;

nVec=params.bundle.nVec;

updateFreq=params.bundle.updateFreq;





breakpoints=params.patternInfo.breakpoints;

x=params.arrayInfo.x;
xSti=params.patternInfo.xStim;



els=[];
for n=1:length(templates)
    
    spikes{n}=NaN*zeros(maxCond,size(TracesAll,2));
    [a b]=sort(max(abs(templates{n}')),'descend');
    ind2=find(a>thresEI);
    if(isempty(ind2))
        ind2=1;
    else
        
    end
    if((length(ind2)==1)&&(b(1)==patternNo))
            ind2=[1 2];
    end
    
    params.neuronInfo.ActiveElectrodes{n}=b(ind2);
    els=union(b(ind2),els);
end

tarray2=setdiff(tarray,0);
params.neuronInfo.ActiveElectrodesAll=els;



elExtra=zeros(size(Art,1),1);
br=intersect(unique(breakpoints{1}(2:end)+1),[1:size(Art,1)]);
elExtra(br)=patternNo;


for i=1:length(elExtra);
    if(elExtra(i)>0)
     
indpat=find(elExtra(i)==els);
   
        indels{i}=setdiff(1:Tmax*length(els),indpat:length(els):Tmax*length(els));
    else
        
        indels{i}=1:Tmax*length(els);
    end
end



KnAll=zeros(length(els)*Tmax,1);
indNeurons=[0 kron([1:length(templates)],ones(1,length(tarray2)))];
indTimes=[0 kron(ones(1,length(templates)),tarray2)];
indTimesArray=[0 kron(ones(1,length(templates)),[1:length(tarray2)])];

for n=1:length(templates)
    for t=1:length(tarray2)
        [ActionPotential]=makeActionPotential(n,tarray2(t),templates,Tmax);
        
        Knn(n,t,:,:)=ActionPotential(:,:);
        Aux(:,t)=reshape(ActionPotential(els,:),Tmax*length(els),1)';
        
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



indStim=find(i<=breakpoints{1}(2:end));
indStim=indStim(1);
indStimOld=indStim;
dLStim=dLSti{indStim};
QStim=QSti{indStim};
QtStim=QtSti{indStim};
KersStim=KersSti{indStim};
xStim=xSti(indStim,:);

krondiag0Stim=1;
krondiag0Stim=kron(krondiag0Stim,dLStim{1});

krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
krondiaginvStim=(exp(xStim(end))*krondiag0Stim*KersStim{2}(i,i)+var0).^(-1);


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
    
    
    
    TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
    
    indTrial=[1:trialI];
    contFindNeurons=1;
    times=zeros(length(templates),trialI);
    while(length(indTrial)>0)
        AA=reshape(TracesResidual(indTrial,els2,1:Tmax),length(indTrial),Tmax*length(els2))'-repmat(AA0,1,length(indTrial));
        
        
        corrs=-2*AA'*KnAll(indels{i},:)+repmat(nansum(KnAll(indels{i},:).^2),length(indTrial),1);
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
        
        TracesResidual(indTrial,:,:)=TracesResidual(indTrial,:,:)-KnnReshaped(idxSubtract,:,:);
        contFindNeurons=contFindNeurons+1;
    end
    
    
    
    
    Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
    
    ArtF(1,ind,:)=FilterArtifactLocal(Kers,Art(1,ind,1:Tmax),[x log(var0)],1,ind,Q,Qt,krondiaginv);
    ArtF(1,patternNo,:)=FilterArtifactLocalStim(KersStim,squeeze(Art(1,patternNo,1:Tmax))',[xStim log(var0)],1,QStim,QtStim,krondiaginvStim);
    
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





xold=x;
KersOld=Kers;
Qold=Q;
dLold=dL;
Qtold=Qt;
krondiag0old=krondiag0;
x01=xold(end);

for i=2:maxCond
    
    
els2=setdiff(els,elExtra(i));
    

i
    if(i>=params.bundle.onsBundle&&params.bundle.onsBundle>1)
        
        if(i==params.bundle.onsBundle&&params.bundle.onsBundle>1)
        
            
            DiagsPre=params.patternInfo.Diags;
            DifsPre=params.patternInfo.Difs;

            DiagsPre{3}=DiagsPre{3}(params.bundle.onsBundle-1);

            DifsPre{3}=DifsPre{3}(params.bundle.onsBundle-1,params.bundle.onsBundle-1);



            f11=@(Art,x)logDetKron(Art(params.bundle.onsBundle-1,ind,:),[xold(1:9)  x log(var0)],DifsPre,[10],[1 1 1],DiagsPre,[3 3 3]);
             g11=@(x)f11(Art,x);
             factor= fminunc(g11,[xold(end)],options)
            x(end)=0;
        
        
    [Res0]=ResidualsElectrodeSimple(Art(1:params.bundle.onsBundle-1,:,:),stimElec,[1:Tmax]);
     %factor=log(nansum(Res0(end,:)')/trace(KersOld{2}))
   
            [Res]=ResidualsElectrodeSimple(Art(params.bundle.onsBundle:end,:,:),stimElec,[1:Tmax]);
           
            %[a bb c]=svd(Res(:,ind));
            %v1=max(a(:,1:nVec)*bb(1:nVec,1:nVec)*c(:,1:nVec)',0);
            [a bb c]=svd(Res(:,ind)-repmat(Res0(end,ind),size(Res,1),1));
            v1=max(a(:,1:nVec)*bb(1:nVec,1:nVec)*c(:,1:nVec)',0)+repmat(Res0(end,ind),size(Res,1),1);
           %v1=max(a(:,1:nVec)*bb(1:nVec,1:nVec)*c(:,1:nVec)',0);
        end
        
        
        Difs{3}=0;
        Diags{3}=1;
        
        DiagsBundle=Diags;
        DiagsBundle{2}=v1(i-params.bundle.onsBundle+1,:)';
        if(mod(i-params.bundle.onsBundle,updateFreq)==0)
           % f11=@(Art,x)logDetKron(Art(i,ind,:),[xold(1:9) x log(var0)],Difs,10,[1 4 1],DiagsBundle,[3 3 3]);
            
            f11=@(Art,x)logDetKron(Art(i,ind,:),[xold(1:3) x xold(4:5) xold(7:9)  0 log(var0)],Difs,[4],[1 5 1],DiagsBundle,[3 3 3],KersOld{2}*exp(factor));
            
            g11=@(x)f11(Art,x);
            
            x11 = fminunc(g11,x01,options)
            x01=[x11];
            x(end)=0;
            
        end
   
        k=2;
       [Ker KerD]=evalKernels(Difs{k},DiagsBundle{k},[x11 xold(4:5)],5,KersOld{2}*exp(factor));
        
        KersNew{k}=Ker;
        
        
        [a b]=eig(KersNew{k});
        Q{k}=a';
        Qt{k}=a;
        dL{k}=diag(b);
        Kers{k}=KersNew{k};
        
        
        krondiag0=kron(dL{1},dL{2});
        
        
        
    end
    


    
    
indStim=find(i<=breakpoints{1}(2:end));
indStim=indStim(1);
breakIni=breakpoints{1}(indStim)+1;
if(~(indStim==indStimOld))

dLStim=dLSti{indStim};
QStim=QSti{indStim};
QtStim=QtSti{indStim};
KersStim=KersSti{indStim};
xStim=xSti(indStim,:);

krondiag0Stim=1;
krondiag0Stim=kron(krondiag0Stim,dLStim{1});

indStimOld=indStim;
end



    krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
    krondiaginvStim=(exp(xStim(end))*krondiag0Stim+var0).^(-1);

    
    krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
    krondiaginvStim=(exp(xStim(end))*krondiag0Stim+var0).^(-1);

    trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
    params.patternInfo.nTrials(i)=trialI;
    
[Apred1]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,i,ArtF(:,ind,:),x,var0);   
Apred(:,ind,:)=Apred1;
if(size(ArtF,1)>=breakIni)
[Apred2]=ExtrapolateArtifactCondStim(KersStim,QStim,QtStim,dLStim,i-breakIni+1,squeeze(ArtF(breakIni:end,patternNo,:)),xStim,var0);   
else
    Apred2=ArtF(end,patternNo,:);
end
Apred(:,patternNo,:)=Apred2;

    flag=1;
    cont=1;
    while(flag==1&&cont<=maxIter)
        
        
        clear times
        
    AA0=reshape(Apred(:,els2,:),Tmax*length(els2),1);
    
    
    
    TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
    
    indTrial=[1:trialI];
    contFindNeurons=1;
    times=zeros(length(templates),trialI);
    
    while(length(indTrial)>0)
        AA=reshape(TracesResidual(indTrial,els2,1:Tmax),length(indTrial),Tmax*length(els2))'-repmat(AA0,1,length(indTrial));
        
        
        corrs=-2*AA'*KnAll(indels{i},:)+repmat(nansum(KnAll(indels{i},:).^2),length(indTrial),1);
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
        
        TracesResidual(indTrial,:,:)=TracesResidual(indTrial,:,:)-KnnReshaped(idxSubtract,:,:);
        contFindNeurons=contFindNeurons+1;
    end
    
    
    
        Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
        
        ArtF(i,ind,:)=FilterArtifactLocal(Kers,Art(1:i,ind,:),[x log(var0)],i,ind,Q,Qt,krondiaginv);
        if(breakIni<i)
        ArtF(i,patternNo,:)=FilterArtifactLocalStim(KersStim,squeeze(Art(breakIni:i,patternNo,:)),[xStim log(var0)],i-breakIni+1,QStim,QtStim,krondiaginvStim);
        else
        ArtF(i,patternNo,:)=FilterArtifactLocalStim(KersStim,reshape(squeeze(Art(breakIni:i,patternNo,:)),1,Tmax),[xStim log(var0)],i-breakIni+1,QStim,QtStim,krondiaginvStim);
        end   
        
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
